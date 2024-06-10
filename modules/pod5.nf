process CONVERT_FAST5_TO_POD5 {
    label 'cpu_8'
    label 'mem_8'
    label 'time_30m'

    container 'quay.io/sangerpathogens/pod5:0.3.6'
    
    input:
    path(fast5s)

    output:
    path("converted.pod5"), emit: pod5_ch

    script:
    """
    pod5 convert fast5 --output converted.pod5 -t ${task.cpus} ${fast5s}
    """
}

process MERGE_POD5 {
    label 'cpu_8'
    label 'mem_8'
    label 'time_30m'

    container 'quay.io/sangerpathogens/pod5:0.3.6'
    
    input:
    path(pod5s)

    output:
    path("merged.pod5"), emit: pod5_ch

    script:
    """
    pod5 merge --output merged.pod5 -t ${task.cpus} ${pod5s}
    """
}