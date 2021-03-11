# Commands to generate test data through the pipeline

## Null models

### Gaussian family

```bash
null_model.py --cluster_file=cluster_cfg.json testdata/null_model.config
cp data/test_null_model_invnorm.RData testdata/null_model.RData
cp data/test_null_model_invnorm_reportonly.RData testdata/null_model_reportonly.RData
```

### Binary

```bash
null_model.py --cluster_file=cluster_cfg.json testdata/null_model_binary.config
cp data/test_null_model.RData testdata/null_model_binary.RData
cp data/test_null_model_reportonly.RData testdata/null_model_binary_reportonly.RData
```

### Unrelated only

```bash
null_model.py --cluster_file=cluster_cfg.json testdata/null_model_unrel.config
cp data/test_null_model_invnorm.RData testdata/null_model_unrel.RData
cp data/test_null_model_invnorm_reportonly.RData testdata/null_model_unrel_reportonly.RData
```
