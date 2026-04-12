# External Dependencies

This directory is intentionally empty in the standalone `TreeCint` repo checkout.
Initialize dependencies as git submodules:

```bash
git submodule update --init --recursive
```

If submodules are not yet added in this clone, add them first:

```bash
git submodule add https://github.com/sunqm/libcint external/libcint
git submodule add https://github.com/fastalgorithms/libid external/libid
git submodule add https://github.com/flatironinstitute/dmk.git external/dmk
git submodule add https://github.com/haiszhu/treefun.git external/treefun
```
