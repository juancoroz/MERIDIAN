# Third-Party Licenses

This document lists the third-party libraries that the MERIDIAN
project depends on. These libraries are **not** distributed with
this source tree; users download them per
`rocket_6dof/BUILD_LINUX.md` or `rocket_6dof/BUILD_WINDOWS.md` and
place them under `rocket_6dof/third_party/`. They are included here
purely for transparency so engineers who populate `third_party/`
know what licenses they are picking up.

The MERIDIAN project's own license (covering the rocket simulator,
the OSK-style kernel, and the example programs) is in `LICENSE` at
the top of this repository.

---

## Optional dependencies (required for `make mission_gui`)

The command-line binaries (`mission`, `mission_mc`, `mission_sobol`,
the regression tests) build with no third-party dependencies.
The graphical interface `mission_gui` requires the three libraries
below.

### Dear ImGui

- **License:** MIT (SPDX: `MIT`)
- **Copyright:** Omar Cornut
- **Upstream:** https://github.com/ocornut/imgui
- **Tested version:** v1.91.6
- **Expected location:** `rocket_6dof/third_party/imgui/`
- **Required by:** `mission_gui`

```
The MIT License (MIT)

Copyright (c) 2014-2024 Omar Cornut

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
```

---

### ImPlot

- **License:** MIT (SPDX: `MIT`)
- **Copyright:** Evan Pezent
- **Upstream:** https://github.com/epezent/implot
- **Tested version:** v0.16
- **Expected location:** `rocket_6dof/third_party/implot/`
- **Required by:** `mission_gui` (plotting panels)

```
MIT License

Copyright (c) 2020 Evan Pezent

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

---

### stb_image_write

- **License:** Dual-licensed: MIT *or* Public Domain (user's choice)
- **Copyright:** Sean Barrett
- **Upstream:** https://github.com/nothings/stb
- **Tested version:** v1.16 (single header: `stb_image_write.h`)
- **Expected location:** `rocket_6dof/third_party/stb/stb_image_write.h`
- **Required by:** `mission_gui` (screenshot export)

Sean Barrett's `stb` libraries are released under a dual license; pick
whichever applies to your use case.

```
This software is available under 2 licenses -- choose whichever you prefer.
------------------------------------------------------------------------------
ALTERNATIVE A - MIT License
Copyright (c) 2017 Sean Barrett
Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
------------------------------------------------------------------------------
ALTERNATIVE B - Public Domain (www.unlicense.org)
This is free and unencumbered software released into the public domain.
Anyone is free to copy, modify, publish, use, compile, sell, or distribute this
software, either in source code form or as a compiled binary, for any purpose,
commercial or non-commercial, and by any means.
In jurisdictions that recognize copyright laws, the author or authors of this
software dedicate any copyright interest in the software to the public domain.
We make this dedication for the benefit of the public at large and to the
detriment of our heirs and successors. We intend this dedication to be an overt
act of relinquishment in perpetuity of all present and future rights to this
software under copyright law.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
```

---

## Reference works (not redistributed)

The MERIDIAN simulator's physics and subsystem models follow the
modeling approach of:

- Peter H. Zipfel, *Modeling and Simulation of Aerospace Vehicle
  Dynamics*, Third Edition, AIAA, 2014.
- The accompanying CADAC++ source code distributed with that work.

These are copyrighted by their respective owners. MERIDIAN's
implementation is an independent reimplementation based on the
published modeling approach; no Zipfel or CADAC++ source code is
included in this repository. See the project README for details.

The OSK-style kernel in `osk/` is independently reimplemented from
the design described in Wilkerson and Smith, *CMD User's Guide*
(Technical Report TR-AMR-SG-05-12, DTIC accession ADA433836, 2005),
an open-source U.S. Government technical report. No code from the
original (non-public) OSK by Ray Sells is included.

---

## Updating this document

When upgrading a vendored library:

1. Verify the upstream license has not changed (most projects use
   the same license across versions, but always check the
   `LICENSE` file in the new release).
2. Update the **Tested version** line above.
3. If the upstream `LICENSE` file's copyright year range has
   advanced, update the inline license text to match.
