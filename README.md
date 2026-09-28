# 10-3rd-party/OCCT-7.8.1/ — mathstudio vendored OpenCASCADE 7.8.1

> OCCT 7.8.1 is now vendored as a mathstudio 3rd-party (this directory).
> Build it once with the script below; the rest of mathstudio consumes it via
> `simplecax-external-kernel::external_kernel_vxloft14` (which links
> `20-self-maintained/03-products/02-od/vendor/occt-config.cmake`).

This directory is a **git submodule** pointing at
[`plasmayang/OCCT`](https://github.com/plasmayang/OCCT), branch
`develop` (pinned to the upstream `V7_8_1` tag — `bd2a789f15235755ce4d1a3b07379a2e062fdc2e`).

The `-patched` suffix is intentionally NOT used here because OCCT fork
contains **no mathstudio-specific patches** — it is a plain upstream fork
kept under the `plasmayang/` org for consistency with
`Catch2-2.11.3-patched/`. If you do apply mathstudio-specific patches to
OCCT, please rename this directory to `OCCT-7.8.1-patched/` and update
the references in `occt-config.cmake` + `.gitmodules`.

---

## Build (Windows / MSVC v142)

```powershell
# From the mathstudio repo root:
10-3rd-party\OCCT-7.8.1\build_occt.bat
```

This will:
1. Pin MSVC v142 toolset (cl 19.29.x) via `vswhere` + `vcvars64.bat`.
2. Configure with `cmake -G Ninja` (or `NMake Makefiles` if Ninja is absent).
4. Build only the modules required for STEP support — FoundationClasses,
   ModelingData, ModelingAlgorithms, DataExchange. Visualization / Draw /
   ApplicationFramework are OFF (headless).
5. Install to `10-3rd-party/OCCT-7.8.1/install/`.
   - Headers: `install/include/opencascade/*.hxx`
   - Libs:    `install/lib/*.lib` (MSVC) / `*.a` (MinGW)

Idempotent: skips configure if `install\include\opencascade\Geom_BSplineCurve.hxx`
already exists. Delete the `install/` tree to force a rebuild.

Approximate build time: **~30 minutes** on a `-j 4` machine. The build log
goes to stdout; pipe to a file if you want to inspect it after.

## Build (Linux)

Use the same OCCT submodule but with the standard cmake workflow:

```bash
cd 10-3rd-party/OCCT-7.8.1
cmake -S . -B build -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
      -DCMAKE_INSTALL_PREFIX="$(pwd)/install" \
      -DBUILD_LIBRARY_TYPE=Static \
      -DBUILD_SHARED_LIBS=OFF \
      -DBUILD_EXAMPLES=OFF \
      -DBUILD_TESTING=OFF \
      -DUSE_TCL=OFF -DUSE_TK=OFF -DUSE_FREETYPE=OFF \
      -DUSE_OPENGL=OFF -DUSE_GLES2=OFF \
      -DBUILD_MODULE_Visualization=OFF \
      -DBUILD_MODULE_ApplicationFramework=OFF \
      -DBUILD_MODULE_Draw=OFF \
      -DBUILD_MODULE_DETools=OFF \
      -DBUILD_MODULE_FoundationClasses=ON \
      -DBUILD_MODULE_ModelingData=ON \
      -DBUILD_MODULE_ModelingAlgorithms=ON \
      -DBUILD_MODULE_DataExchange=ON
cmake --build build -j 4
cmake --install build
```

## Consumed by mathstudio

After building, the `20-self-maintained/03-products/02-od/vendor/occt-config.cmake`
file (included unconditionally by `02-od/CMakeLists.txt`) discovers OCCT in
this directory's `install/` subtree and exposes the libraries as
`OpenCASCADE::*` imported targets. `simplecax-external-kernel::external_kernel_vxloft14`
(the VxLoft14 wrapper used by `e2e-gallery/od/` and `e2e-gallery/native/`)
then links against those.

## Pinning strategy

The submodule is pinned to:
- URL: `https://github.com/plasmayang/OCCT.git`
- Branch: `develop` (mathstud matrix-tracked; **NOT** a synced default
  branch from upstream)
- Commit: `bd2a789f15235755ce4d1a3b07379a2e062fdc2e` (upstream V7_8_1)

If you need to bump to a newer OCCT release:
1. Create a new branch on the fork: `git push origin V<X>_<Y>_<Z>:<new-branch>`
2. Update `.gitmodules` branch + the parent submodule commit pointer.
3. Update `occt-config.cmake`'s expected version comment if the API drifts.
4. Update this README.