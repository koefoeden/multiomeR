rlang::list2(
  tarchetypes::tar_file(
    name = CollecTRI_human_network_csv,
    description = "Download and checksum the published human CollecTRI signed TF-target network [part_of_graph:differential_analyses]",
    command = download_checksum_verified_file(
      url = collectri_source()$url,
      expected_sha256 = collectri_source()$sha256,
      destfile = file.path(targets::tar_config_get("store"), "files", "CollecTRI", "CollecTRI_human_2023.csv")
    )
  ),
  targets::tar_target(
    name = CollecTRI_human_network_tibble,
    description = "Validate and reformat the checksum-pinned human CollecTRI signed TF-target network [part_of_graph:differential_analyses]",
    command = read_CollecTRI_human_network(CollecTRI_human_network_csv)
  )
)
