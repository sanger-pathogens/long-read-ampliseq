def mark_read_duplicates_in_summary(sequencing_summary, outputFilePath, keep_or_remove){
    def filePath = sequencing_summary.toString()

    // Read the TSV file
    def tsvFile = new File(filePath)

    //make a set to store what we have seen
    def column2Set = new HashSet<String>()
    def duplicates=[]

    // Iterate over each line in the file
    tsvFile.eachLine { line ->
        // Split the line by tabs
        def columns = line.split('\t')
    
        // Ensure the line has at least two columns
        if (columns.size() >= 2) {
            def value = columns[1]
            if (column2Set.contains(value)) {
                duplicates << value
            } else {
                column2Set.add(value)
            }
        }
    }

    //for samtools view if the file starts with ^ it removes all entrys from the file instead of keeping them
    def finalPath = "${outputFilePath}/duplicates_${workflow.runName}.txt"

    def outputFile = new File(finalPath)

    duplicates.each { value ->
        outputFile.append(value + '\n')
    }

    return finalPath
}

process SUMMARY_DUPLICATES { 
    label 'cpu_1'
    label 'mem_4'
    label 'time_1'
    
    input:
    val(summary)
    val(mode)

    output:
    path("*.txt"), emit: summary_channel

    exec:
    mark_read_duplicates_in_summary(summary, "${task.workDir}", mode)
}