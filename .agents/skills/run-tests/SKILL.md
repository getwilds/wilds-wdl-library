---
name: run-tests
description: Run the testrun.wdl for a WILDS module or pipeline using sprocket or miniwdl
argument-hint: <ww-name> [executor]
allowed-tools: Bash(make *), Read, Glob, Grep
---

# Run Tests for a WILDS Module or Pipeline

Execute the `testrun.wdl` for `$0` and report results.

## Steps

### 1. Locate the Test Workflow

Find `testrun.wdl` in either:
- `modules/$0/testrun.wdl`
- `pipelines/$0/testrun.wdl`

### 2. Run the Test

If a second argument (`$1`) is provided, use that executor. Otherwise default to sprocket.

**Sprocket:**
```bash
make run_sprocket NAME=$0
```

**miniWDL:**
```bash
make run_miniwdl NAME=$0
```

### 3. Report Results

- If the test passes, report success with a brief summary of outputs
- If the test fails, analyze the error output, suggest fixes, and offer to apply them
