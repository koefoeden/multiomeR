# CollecTRI human regulons

The public pipeline downloads the published human CollecTRI signed
transcription-factor target network at runtime from the OmniPath rescue
archive:

```text
https://rescued.omnipathdb.org/CollecTRI.csv
```

SHA-256:

```text
86c90b30f2cc75c189da1f0a8c353d1547287cd656a9fac1c678634285bcb4e0
```

The download target rejects content that does not match this checksum. The
validated snapshot contains 43,536 unique signed source-target interactions
for 1,189 transcription-factor or TF-complex sources. AP1 and NFKB are retained
as complex regulons. CollecTRI is described in Müller-Dott et al., *Nucleic
Acids Research* 2023, https://doi.org/10.1093/nar/gkad841.

## Redistribution boundary

OmniPath classifies CollecTRI as a composite resource whose constituent
databases retain their own licenses. The public repository therefore contains
the checksum-pinned download and validation code, not a redistributed copy of
the network. Users remain responsible for complying with the applicable
constituent licenses.
