# Source and license notices

The comparison toolbox upstream README specifies GPL v3:
https://github.com/Baoshan-Song/KFV-FGO-Comparison
The Python port preserves this attribution; see LICENSE for GPL v3 text.

Real GNSS/IMU parsing and propagation were adapted from ZzPolyU/KF-GNSS-INS.
MatRTKLIB's wrapper conventions inform the Python GNSS adapter; original notices
for both MIT-licensed sources are retained below. The small WGS84 coordinate
conversion follows RTKLIB/MALIB; its notice is included as well.

pyRTKLIB is an installed dependency, not a bundled native binary in this source
distribution: https://github.com/IPNL-POLYU/pyrtklib. Its own distribution contains
its binding and RTKLIB notices. NumPy, SciPy and Matplotlib are installed separately.

Example GNSS/IMU data are copied unchanged from the existing user-provided
urban navigation recording in the MATLAB comparison project. No new license or
ownership claim is made over the recording or its ground truth.

## KF-GNSS-INS

MIT License

Copyright (c) 2025 Zhan Zhi

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

## MatRTKLIB

MIT License

Copyright (c) 2023 Taro Suzuki

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

## RTKLIB / MALIB


The MALIB software package is distributed under the following BSD 2-clause
license. Users are permitted to develop, produce or sell their own non-
commercial or commercial products utilizing, linking or including MALIB as long
as they comply with the license.

--------------------------------------------------------------------------------

         Copyright (c) 2007-2023, T. Takasu, All rights reserved.
         Copyright (c) 2023-2024, Japan Aerospace Exploration Agency. All Rights Reserved.
         Copyright (C) 2023-2024, TOSHIBA ELECTRONIC TECHNOLOGIES CORPORATION. All Rights Reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list
of conditions and the following disclaimer. Redistributions in binary form must
reproduce the above copyright notice, this list of conditions and the following
disclaimer in the documentation and/or other materials provided with the
distribution.

The software package includes some companion executive binaries or shared
libraries necessary to execute APs on Windows. These licenses succeed to the
original ones of these software.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

--------------------------------------------------------------------------------

Notes:
Previous versions of RTKLIB until ver. 2.4.1 had been distributed under GPLv3
license.
