rlang::list2(
  tarchetypes::tar_file(
    name = CollecTRI_human_network_csv,
    description = "Download and checksum the published human CollecTRI signed TF-target network [part_of_graph:differential_analyses]",
    command = download_CollecTRI_human_network(
      network_url = collectri_source()$url,
      expected_sha256 = collectri_source()$sha256
    )
  ),
  targets::tar_target(
    name = CollecTRI_human_network_tibble,
    description = "Validate and reformat the checksum-pinned human CollecTRI signed TF-target network [part_of_graph:differential_analyses]",
    command = read_CollecTRI_human_network(CollecTRI_human_network_csv)
  )
)
