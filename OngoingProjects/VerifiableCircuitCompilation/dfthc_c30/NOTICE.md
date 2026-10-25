# Attribution

This standalone loader adapts `gauge_qrom.py` and the required parts of
`banked_lab.py` from Dmitrii Gluzdov's Apache-2.0 public notebook,
[version 8](https://www.kaggle.com/code/dmitriigluzdov/quantum-circuit-compilation-from-source).
The notebook's kernel ID is 140416287; its source and module hashes are pinned
in the sibling [release audit](../dfthc_feasibility/public-release-provenance.json)
and [loader audit](../dfthc_feasibility/released-loader.json).
The standalone port removes the notebook modules and adds validation,
reproduction scripts, WL identities and full-verifier evidence.

The banked QROM construction credits Kozykorpesh / Sumit Chongder:
[Kozykappa source](https://www.kaggle.com/code/kozykappa/quantum-ai-compact-arithmetic-banked-qrom?scriptVersionId=361433617),
[Sumit Chongder source](https://www.kaggle.com/code/sumitchongder/notebook2f196e646b?scriptVersionId=361357149),
as read through Kirill Tym's public quantum-circuits-qrom-tower42 v10 source.
The signed gauge, four-address interpolation, selective AND pebbling and
split-ANF measurement cleanup are attributed public work.

The fixed factorization comes from the Apache-2.0
[FeMoCo public factorization, version 2](https://www.kaggle.com/datasets/dmitriigluzdov/vqcc-femoco-public-factorization/versions/2).
Its SHA256 is
`9afd0ca5e90dd1740fe713e4ac1c26736184ea9bc4d4b57a6eeefafd608ee193`.
Downloaded factors and generated circuits are excluded from public files.

Source code in this folder is provided under the accompanying [Apache-2.0 license](LICENSE).
