using CellularPotts, Literate, Random 


Random.seed!(314159)

#Function to replace name of gif with actual gif
function gifReplace(content, root, fileName)

    gifFile = replace(fileName, ".jl" => ".gif")

    #Gotta love Windows
    partialPath = replace(root, "./"=>"", "\\"=>"/")
    newURLPath = join([partialPath, gifFile], "/")

    oldStr = "````\n\"$(gifFile)\"\n````"

    newStr = 
    """
    ```@raw html
    <p style="text-align:center;">
        <img title="" src="https://github.com/RobertGregg/CellularPotts.jl/blob/master/$(newURLPath)?raw=true" width="445">
    </p>
    ```
    """

    return replace(content, oldStr => newStr)
end

#Loop through all the examples, execute them, save markdown file

examplesToUpdate = ["OverHere.jl"]

for (root, dirs, files) in walkdir("./docs/src/ExampleGallery")
    for file in files
        if endswith(file,".jl") && file ∈ examplesToUpdate
            Literate.markdown(joinpath(root, file), root; execute=true, postprocess=str->gifReplace(str,root,file))
        end
    end
end