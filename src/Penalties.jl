#Note, the Types definitions are in Core.jl (e.g. struct AdhesionPenalty ... end)

####################################################
# Adhesion
####################################################

function addPenalty!(cpm::CellPotts, AP::AdhesionPenalty)
    return addPenalty!(cpm, AP, cpm.step.sourceCellID) - addPenalty!(cpm, AP, cpm.step.targetCellID) 
end

function addPenalty!(cpm::CellPotts, AP::AdhesionPenalty, σᵢ::T) where T<:Integer
    #Initialize the penality
    adhesion = zero(T)

    τᵢ = cpm.currentState.typeIDs[σᵢ]

    for neighbor in cpm.step.targetNeighborNodes
        #Given a node index, get the cellID
        σⱼ = cpm.space.nodeIDs[neighbor]

        #Convert the cellID to cellType
        τⱼ = cpm.currentState.typeIDs[σⱼ]

        #Adhesion is increased if adjacent cells are different types
        adhesion += AP.J[τᵢ, τⱼ] * (1-δ(σᵢ, σⱼ))
    end

    return adhesion
end

####################################################
# Volume
####################################################

function addPenalty!(cpm::CellPotts, VP::VolumePenalty)
    σᵢ = cpm.step.sourceCellID
    σⱼ = cpm.step.targetCellID
    sourceVolume = addPenalty!(cpm, VP, σᵢ) + addPenalty!(cpm, VP, σⱼ)
   
    #Change the volumes and recalculate penalty
    cpm.currentState.volumes[σᵢ] += 1
    cpm.currentState.volumes[σⱼ] -= 1

    targetVolume = addPenalty!(cpm, VP, σᵢ) + addPenalty!(cpm, VP, σⱼ)

    #Reset the volumes
    cpm.currentState.volumes[σᵢ] -= 1
    cpm.currentState.volumes[σⱼ] += 1

    return targetVolume - sourceVolume
end


function addPenalty!(cpm::CellPotts, VP::VolumePenalty, σ::T) where T<:Integer
    volume = cpm.currentState.volumes[σ]
    desiredVolume = cpm.currentState.desiredVolumes[σ]
    τⱼ = cpm.currentState.typeIDs[σ]

    return VP.λᵥ[τⱼ] * (volume - desiredVolume)^2
end


####################################################
# Perimeter
####################################################

function addPenalty!(cpm::CellPotts, PP::PerimeterPenalty)
    node = cpm.step.targetNode
    σᵢ = cpm.step.sourceCellID
    σⱼ = cpm.step.targetCellID
    sourcePerimeter = addPenalty!(cpm, PP, σᵢ) + addPenalty!(cpm, PP, σⱼ)

    #Unlike volumes which change by ±1, perimeter is more complicated
    PP.Δpᵢ = -perimeterLocal(cpm.space, node, σᵢ)
    PP.Δpⱼ = -perimeterLocal(cpm.space, node, σⱼ)

    #Copy the source node into the target node
    cpm.space.nodeIDs[node] = σᵢ

    PP.Δpᵢ += perimeterLocal(cpm.space, node, σᵢ)
    PP.Δpⱼ += perimeterLocal(cpm.space, node, σⱼ)

   
    #Change the perimeters and recalculate penalty
    cpm.currentState.perimeters[σᵢ] += PP.Δpᵢ
    cpm.currentState.perimeters[σⱼ] -= PP.Δpⱼ

    targetPerimeter = addPenalty!(cpm, PP, σᵢ) + addPenalty!(cpm, PP, σⱼ)

    #Reset the perimeters (and space)
    cpm.currentState.perimeters[σᵢ] -= PP.Δpᵢ
    cpm.currentState.perimeters[σⱼ] += PP.Δpⱼ

    cpm.space.nodeIDs[node] = σⱼ

    return targetPerimeter - sourcePerimeter
end

function addPenalty!(cpm::CellPotts, PP::PerimeterPenalty, σ::T) where T<:Integer

    perimeter = cpm.currentState.perimeters[σ]
    desiredPerimeter = cpm.currentState.desiredPerimeters[σ]
    τⱼ = cpm.currentState.typeIDs[σ]

    return PP.λₚ[τⱼ] * (perimeter - desiredPerimeter)^2
end


function perimeterLocal(space::CellSpace, n₀::T, σ::T) where T<: Integer
    
    perimeter = zero(T)

    #Loop through the neighbors and calculate each neighbor perimeter penality
    for n₁ in neighbors(space,n₀)
        #Only add perimeter penality for the given cellID
        if space.nodeIDs[n₁] == σ
            #To calculate each neighbor perimeter penality, loop through it's neighbors
            for n₂ in neighbors(space,n₁)
                #Penalize Perimeter if the pixels do not match
                if space.nodeIDs[n₂] ≠ space.nodeIDs[n₁]
                    perimeter += 1
                end
            end
        end
    end

    #The source/target node only contributes if the current space node matches
    if space.nodeIDs[n₀] == σ
        for n₁ in neighbors(space,n₀)
            if space.nodeIDs[n₁] ≠ space.nodeIDs[n₀]
                perimeter += 1
            end
        end
    end

    return perimeter
end

####################################################
# Migration
####################################################

#=
The orginal method to calculate migration penality uses geometric mean. This is problematic because it introduces floating point calculations. Here we instead calculate an integer rounded arithmetic mean and map any neighorbood containing zeros to zero. This mimics the desired behavior of the geometric mean where neighborhoods “with holes” (i.e., lattice sites with activity value zero) are ignored.

This change might become moot when chemotaxis is introduced. Gradient fields with unavoidably introduce floating point calculations.
=#

function addPenalty!(cpm::CellPotts, MP::MigrationPenalty)

    return MP.λ * (addPenalty!(cpm, MP, cpm.step.targetNode) - addPenalty!(cpm, MP, cpm.step.sourceNode))
end


function addPenalty!(cpm::CellPotts, MP::MigrationPenalty, σ::T) where T<:Integer

    memoryNeighbors = view(MP.nodeMemory, neighbors(cpm.space,σ))
    numNeighbors = length(memoryNeighbors)
    return any(iszero,memoryNeighbors) ? 0 : sum(memoryNeighbors) ÷ numNeighbors
end

####################################################
# Chemotaxis
####################################################