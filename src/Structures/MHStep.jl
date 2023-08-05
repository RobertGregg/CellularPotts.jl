####################################################
# Variables for Markov Step 
####################################################

mutable struct CandidateNode{T<:Integer}
    node::T                 #Index of node choosen
    neighbors::Vector{T}    #Indicies for the neighboring nodes
    id::T                   #ID of choosen cell
    type::T                 #Type of the choosen cell
end

CandidateNode() = CandidateNode(0,[0],0,0)

mutable struct MHStep{T<:Integer}
    source::CandidateNode{T}
    target::CandidateNode{T}
    counter::T                #Counts the number of ModelSteps performed
    success::Bool             #Tracks if the MHStep was successful
end 

MHStep() = MHStep(CandidateNode(), CandidateNode(), 0, false)