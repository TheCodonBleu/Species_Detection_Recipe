version 1.0

workflow SpeciesDetection {
    meta {
      description: "Determine your species from the reads sequences."
      author: "David Maimoun (The Codon Bleu)"
  }
    input {
        String basespace_collection_id 
        String api_server
        String access_token
        String? sample_prefix

    }

    call GetReadsList {
        input:
            basespace_collection_id = basespace_collection_id,
            access_token = access_token,
            api_server = api_server,
            sample_prefix = select_first([sample_prefix, ""])

    }

    scatter(sample_name in GetReadsList.samples_name) {
        call FetchReads {
          input:
            basespace_sample_name = sample_name,
            basespace_collection_id = basespace_collection_id,
            api_server = api_server,
            access_token = access_token
        }

        call Detect_Species {
            input:
                read1 = FetchReads.read1,
                read2 = FetchReads.read2,
                sample_id = sample_name
        }
        
    }

    call MergeReports {
        input:
            species_detected_list = Detect_Species.sample_detected
    }

    output {
        File reads_list = GetReadsList.reads_list
        Array[String] samples_name = GetReadsList.samples_name
        File species_detected = MergeReports.species_detected
    
    }

}
  
task GetReadsList {
    input {  
        String basespace_collection_id
        String api_server 
        String access_token
        String? sample_prefix
        String docker = "us-docker.pkg.dev/general-theiagen/theiagen/basespace_cli:1.2.1"

    }

    command <<<     
        bs project content --name ~{basespace_collection_id} \
            --api-server=~{api_server} \
            --access-token=~{access_token} \
            --retry > reads_list.txt

        # Create a a samples name list from the reads_list
        # -> Fetch the name form the R1 read, and remove if got a 'Undetermined' read
        if [ -z ~{sample_prefix} ]; then
            grep -o "[A-Za-z0-9_-]*_S[0-9]*_L[0-9]*_R1_[0-9]*\.fastq\.gz" reads_list.txt \
            | sed 's/_S[0-9]*_L[0-9]*_R1_.*\.fastq\.gz//' \
            | grep -v "^Undetermined$" \
            > samples_name.txt
        else
            grep -o "~{sample_prefix}[A-Za-z0-9_-]*_S[0-9]*_L[0-9]*_R1_[0-9]*\.fastq\.gz" reads_list.txt \
            | sed 's/_S[0-9]*_L[0-9]*_R1_.*\.fastq\.gz//' \
            > samples_name.txt
        fi
    >>>

    output {
        File reads_list = "reads_list.txt"
        Array[String] samples_name = read_lines("samples_name.txt")
    }

    runtime {
        docker: docker
        preemptible: 1
  }
}

task FetchReads {
    input {
        String basespace_sample_name
        String? basespace_sample_id   
        String basespace_collection_id
        String api_server 
        String access_token
     
        String docker = "us-docker.pkg.dev/general-theiagen/theiagen/basespace_cli:1.2.1"
    }

    command <<<
        # Trim basespace_collection_id (run_name) because spaces cause problem
        run_name_trimmed=$(echo "~{basespace_collection_id}" | awk '{$1=$1;print}')

        # set basespace name and id variables
        if [[ ! -z "~{basespace_sample_id}" ]]; then
            sample_identifier="~{basespace_sample_name}"
            dataset_name="~{basespace_sample_id}"
        else
            sample_identifier="~{basespace_sample_name}"
            dataset_name="~{basespace_sample_name}"
        fi
    
        # print all relevant input variables to stdout
        echo -e "sample_identifier: ${sample_identifier}\ndataset_name: ${dataset_name}\nbasespace_collection_id: $run_name_trimmed"
        
        #Set BaseSpace comand prefix
        bs_command="bs --api-server=~{api_server} --access-token=~{access_token}"
        echo "bs_command: ${bs_command}"

        #Grab BaseSpace Run_ID from given BaseSpace Run Name
        run_id=$(${bs_command} list run --retry | grep "$run_name_trimmed" | awk -F "|" '{ print $3 }' | awk '{$1=$1;print}' )
        echo "run_id: ${run_id}" 
        
        if [[ ! -z "${run_id}" ]]; then 
            #Grab BaseSpace Dataset ID from dataset lists within given run 
            dataset_id_array=($(${bs_command} list dataset --retry --input-run=${run_id} | grep "${dataset_name}" | awk -F "|" '{ print $3 }' )) 
            echo "dataset_id: ${dataset_id_array[*]}"
        
        else 
            #Try Grabbing BaseSpace Dataset ID from project name
            echo "Could not locate a run_id via Basespace runs, attempting to search Basespace projects now..."
            
            project_id=$(${bs_command} list project --retry | grep "$run_name_trimmed" | awk -F "|" '{ print $3 }' | awk '{$1=$1;print}' )
            
            echo "project_id: ${project_id}" 

            if [[ ! -z "${project_id}" ]]; then 
                echo "project_id identified via Basespace, now searching for dataset_id within project_id ${project_id}..."
                dataset_id_array=($(${bs_command} list dataset --retry --project-id=${project_id} | grep "${dataset_name}" | awk -F "|" '{ print $3 }' ))  
                echo "dataset_id: ${dataset_id_array[*]}"
            else       
                echo "No run or project id found associated with input basespace_collection_id: $run_name_trimmed" >&2
                exit 1
            fi      
        fi

        #Download reads by dataset ID
        for index in ${!dataset_id_array[@]}; do
            dataset_id=${dataset_id_array[$index]}
            mkdir ./dataset_${dataset_id} && cd ./dataset_${dataset_id}
        
            echo "dataset download: ${bs_command} download dataset -i ${dataset_id} -o . --retry"
            ${bs_command} download dataset --retry -i ${dataset_id} -o . --retry && cd ..
            echo -e "downloaded data: \n $(ls ./dataset_*/*)"
        done

        # rename FASTQ files to add back in underscores that Illumina/Basespace changed into hyphens
        echo "Concatenating and renaming FASTQ files to add back underscores in basespace_sample_name"
        # setting a new bash variable to use for renaming during concatenation of FASTQs
        # SAMPLENAME_HYPHEN_INSTEAD_OF_UNDERSCORES=$(echo $sample_identifier | sed 's|_|-|g' | sed 's|\.|-|g')

        # echo $SAMPLENAME_HYPHEN_INSTEAD_OF_UNDERSCORES > 'sample_id.txt'

        echo $sample_identifier > 'sample_id.txt'
 

        for fwd_read in ./dataset_*/${sample_identifier}_*R1_*.fastq.gz; do
            if [[ -s $fwd_read ]]; then
                read1_name=$(basename "$fwd_read")

                echo ${read1_name} > read1_name.txt
                cat $fwd_read      > fwd.fastq.gz
                
            fi
        done

        for rev_read in ./dataset_*/${sample_identifier}_*R2_*.fastq.gz; do
            if [[ -s $rev_read ]]; then
                read2_name=$(basename "$rev_read")

                echo ${read2_name} > read2_name.txt
                cat $rev_read      > rev.fastq.gz
            fi
        done
    >>>

    output {
        File? read1          = 'fwd.fastq.gz'
        File? read2          = 'rev.fastq.gz'
    }

    runtime {
        docker: docker
        maxRetries: 1
        
  }
}


task Detect_Species {
  input {
    File? read1
    File? read2
    String sample_id
    
    String docker = "bioinfomoh/specie_detection:1"
    Int cpu = 20
  
  }

  command <<<
        mode=""
        compressed=""

        if [ ! -s "~{read1}" ]; then
            echo "~{sample_id} not found" > "~{sample_id}.report"
            echo "~{sample_id},not_found,xxx" > species_detected.csv        
        else
            # Check if paired mode should be used
            if ! [ -z "~{read2}" ]; then
                echo "Reads are paired..."
                mode="--paired"
            fi

            # Determine if reads are compressed
            if [[ "~{read1}" == *.gz ]]; then
                echo "Reads are compressed..."
                compressed="--gzip-compressed"
            fi

            # Run Kraken2
            echo "Running Kraken2..."
            kraken2 $mode $compressed --threads "~{cpu}" --use-names --db /app/db/kraken_db \
                --report "~{sample_id}.report" --paired "~{read1}" "~{read2}" --output -

            echo "~{sample_id}" > ~{sample_id}_sample_detected.csv
            
            declare -A species
            species["NM"]="Neisseria Meningitidis"
            species["NG"]="Neisseria Gonorrhoeae"
            species["HI"]="Haemophilus Influenzae"
            species["SH"]="Salmonella"
            species["SO"]="Salmonella"
            species["LC"]="Listeria monocytogenes"
            species["LF"]="Listeria monocytogenes"
            species["SG"]="Shigella"
            species["CA"]="Campylobacter"
            species["VIB"]="Vibrio"
            species["V"]="Vibrio"
            species["EC"]="Escherichia coli"
            species["SA"]="Staphylococcus aureus"
            species["BP"]="Bordetella pertussis"
            species["SP"]="Streptococcus pneumoniae"
            species["ST"]="Streptococcus pyogenes"
            species["ST"]="Streptococcus agalactiae"
            species["LG"]="Legionella pneumophila"
            species["LW"]="Legionella pneumophila"
            species["CB"]="Corynebacterium diphtheriae"
            species["HI"]="Haemophilus influenzae"
            species["NM"]="Neisseria meningitidis"
            species["M" ]="Neisseria meningitidis" 

            prefix=$(echo "~{sample_id}" | grep -o '^[^0-9]*')

            # Extract detected species from report
            detected=$(awk -F'\t' '$4 == "S" {gsub(/^[ \t]+/, "", $6); print $6; exit}' "~{sample_id}.report")

            # Convert both to lowercase for case-insensitive comparison
            detected_lower=$(echo "$detected" | tr '[:upper:]' '[:lower:]')
            expected_lower=$(echo "${species[$prefix]}" | tr '[:upper:]' '[:lower:]')

            if [[ "$prefix" == "ST" && "$detected" == *"Streptococcus"* ]]; then
                # handle ST samples (can be pyogenes or agalactiae)
                echo "Met ST , and detected Streptococcus"
                echo "~{sample_id},${detected},+" > ~{sample_id}_sample_detected.csv

            elif [[ "$prefix" == "SG" && "$detected" == "Escherichia coli" ]]; then
                # handle Shigella (SG) - Kraken detect it as EC
                echo "~{sample_id},${detected},??" > ~{sample_id}_sample_detected.csv

            elif [[ "$detected_lower" == *"$expected_lower"* ]]; then
                echo "~{sample_id},${detected},+" > ~{sample_id}_sample_detected.csv

            else   
                if [[ -n "${detected}" && "${detected}" != "" ]]; then
                    echo "~{sample_id},${detected},xxx" > ~{sample_id}_sample_detected.csv
                else
                    echo "~{sample_id},,xxx" > ~{sample_id}_sample_detected.csv
                fi
            fi
        fi
    >>>

    output {
        File report = "~{sample_id}.report"
        String sample_detected = read_string("~{sample_id}_sample_detected.csv")
    }

    runtime {
        docker: docker
        cpu: cpu
    }
}

task MergeReports {
    input {
        Array[String] species_detected_list
    }

    command <<<
        echo "Sample,Detected,Match" > species_detected.csv
        echo "~{sep='\n' species_detected_list}" >> species_detected.csv
        
    >>>

    output {
        File species_detected = "species_detected.csv"
    }

    runtime {
        docker: "ubuntu:20.04"

    }
}
