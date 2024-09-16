using Configurations

"""
    Config

Configuration struct for the EnsemblGeneMapping module.

# Fields
- `gene_col::String`: The name of the column in the expression data that contains the gene symbols. Default is "gene_symbol".
- `mart_id::String`: The ID of the BioMart database. Default is "ensembl".
- `mart_dataset::String`: The dataset to use in the BioMart database. Default is "hsapiens_gene_ensembl".
"""
@option struct Config
    gene_col::String = "gene_symbol"
    mart_id::String = "ensembl"
    mart_dataset::String = "hsapiens_gene_ensembl"
end

function load_config(file_path::String)::Config
    if isempty(file_path)
        @debug "No configuration file provided. Using default configuration."
        return Config()
    end

    # Check if the file exist 
    if !isfile(file_path)
        @error "File not found: $file_path"
        exit(1)
    end

    try
        config_dict = TOML.parsefile(file_path)

        return from_dict(Config, config_dict)
    catch e
        if isa(e, SystemError)
            @error "Something went wront while reading $file_path. File exists but could not be read."
        elseif isa(e, TOML.ParserError)
            @error "Invalid TOML format in $file_path."
        else
            rethrow(e)
        end
        exit(1)
    end
end
