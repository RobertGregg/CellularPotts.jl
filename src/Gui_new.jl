function CellGUI(cpm::CellPotts{2})

    ####################################################
    # Simulation Screen
    ####################################################

    #Start by creating a blank figure
    fig = Figure(resolution = (1200, 1200), backgroundcolor = RGBf0(0.98, 0.98, 0.98))

    #Create observables that will change as the simulation progresses
    timestep = Node(1) #will increase by one every step
    maxLength = 1000 #Window history for plot
    energies = Node(fill(cpm.energy, maxLength)) #vector of past cpm energies

    #The energies plot updates overtime and needs dynamic axes limits
    lim = @lift begin
        xmin = max(0, $timestep - maxLength)
        xmax = max(maxLength, $timestep)
        ymin, ymax =  extrema($energies)

        (xmin, xmax, ymin-1, ymax+1)
    end

    #Name the axes on the figure
    axSim = fig[1, 1] = Axis(fig, title = "Simulation")
    axEnergy = fig[1, 2] = Axis(fig, title = "Energy", limits = lim, height = 300,
                                tellheight = false, valign = :top)
    colsize!(fig.layout, 1, Relative(2/3)) #Make cells heatmap plot relatively larger


    ####################################################
    # Figure 1 (heatmap)
    ####################################################

    #The first plot will show a heatmap of the cell simulation
    # heatmap_node is an array that updates when timestep updates
    heatmap_node = @lift begin
        currentTime = $timestep
        MHStep!(cpm)
        cpm.visual
    end

    #Create the heatmap
    heatmap!(axSim,
            heatmap_node,
            show_axis = false,
            colormap = :Purples) #:Greys_3
    tightlimits!.(axSim)
    hidedecorations!.(axSim) #removes axis numbers


    ####################################################
    # Figure 2 (Energy)
    ####################################################

    #On the 2nd axis, plot the energies
    xvals = @lift( $lim[1]+1:$lim[2] )
    lines!(axEnergy, xvals, energies)


    ####################################################
    # Buttons
    ####################################################

    fig[2,1] = buttongrid = GridLayout(tellwidth = false, halign = :left)
    buttonsLabels = ["▶","■"]

    #Loop through and assign button labels, width, and color
    buttongrid[1,1:length(buttonsLabels)] = [Button(
        fig,
        label = lab,
        width = 70,
        buttoncolor = :grey) for lab in buttonsLabels]

    #Set the buttons to individual variables
    playpause, stop = contents(buttongrid)

    #If the play/pause button is clicked, change the label
    on(playpause.clicks) do clicks
        if isodd(clicks)
            playpause.label = "𝅛𝅛"
        else
            playpause.label = "▶"
        end
    end

    display(fig)

    runsim = true
    stop.clicks[] = 0 #for stop button
    while runsim

        #Is the pause button pushed?
        if playpause.label[] == "▶"
            timestep[] += 1
            appendEnergy!(energies,cpm)
            notify(timestep)
            notify(energies)
        end

        #Has the stop button been pushed?
        if  stop.clicks[] == 1
            runsim = false
            GLMakie.destroy!(GLMakie.global_gl_screen()) #close the window
        end

        sleep(eps())
    end

end

#This is suprisingly better than circshift
function appendEnergy!(energies,cpm)
    popfirst!(energies[])
    push!(energies[], cpm.energy)
end