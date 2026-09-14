rlang::list2(
  tarchetypes::tar_file(
    name = CollecTRI_human_network_csv,
    description = "Download and checksum the published human CollecTRI signed TF-target network [part_of_graph:differential_analyses]",
    command = download_CollecTRI_human_network(
      network_url = "https://rescued.omnipathdb.org/CollecTRI.csv",
      expected_sha256 = "86c90b30f2cc75c189da1f0a8c353d1547287cd656a9fac1c678634285bcb4e0"
    )
  ),
  targets::tar_target(
    name = CollecTRI_human_network_tibble,
    description = "Validate and reformat the checksum-pinned human CollecTRI signed TF-target network [part_of_graph:differential_analyses]",
    command = read_CollecTRI_human_network(CollecTRI_human_network_csv)
  )
)
