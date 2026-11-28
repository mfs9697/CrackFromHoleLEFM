# IAM-UA8-F02: rebuild the proposed new Figure 2 (structured crack-tip mesh) from the packed P21 mesh data.
# 1 verification + vector mesh layers (Python, numpy) -> figures_src/iam_f02/IAM_F02_VERIFICATION.json, figures_src/iam_f02/layers/
# 2 the three PGFPlots panels (XeLaTeX, STIX Two) -> figures/iam_f02/{a_global_mesh,b_structured_region,c_symmetry_detail}.pdf
# 3 the two previews at the printed size (layouts A and B) -> figures/iam_f02/preview_layout_{A,B}.pdf
# The mesh data themselves come from scripts/iam_f02_regenerate_mesh.m (MATLAB, pinned generator; see the report). No FEM.
# Usage (repository root): .\paper\ukrainian_review\scripts\build_iam_f02.ps1
$ErrorActionPreference = 'Stop'
$ua = Split-Path -Parent $PSScriptRoot
$src = Join-Path $ua 'figures_src\iam_f02'
foreach ($f in @('layers\iam_f02_a_global_mesh.pdf', 'layers\iam_f02_b_structured_region.pdf', 'layers\iam_f02_c_symmetry_detail.pdf',
                 'layers\iam_f02_constants.tex', 'IAM_F02_VERIFICATION.json')) {
    $p = Join-Path $src $f; if (Test-Path -LiteralPath $p) { Remove-Item -LiteralPath $p }
}
Get-ChildItem -LiteralPath (Join-Path $src 'layers') -Filter *.csv -ErrorAction SilentlyContinue | Remove-Item
& python (Join-Path $PSScriptRoot 'iam_f02_mesh_figure.py')
if ($LASTEXITCODE -ne 0) { throw 'iam_f02_mesh_figure.py failed (verification or layers)' }
& (Join-Path $PSScriptRoot 'build_ua_figures.ps1') -Only iam_f02/a_global_mesh, iam_f02/b_structured_region, iam_f02/c_symmetry_detail
$work = Join-Path $src '.build'; $null = New-Item -ItemType Directory -Path $work -Force
$dest = Join-Path $ua 'figures\iam_f02'
foreach ($n in @('preview_layout_A', 'preview_layout_B')) {
    Push-Location $src
    try {
        & xelatex -interaction=nonstopmode -halt-on-error -file-line-error "-output-directory=$work" "$n.tex" | Out-Null
        if ($LASTEXITCODE -ne 0) { throw "xelatex failed for $n (see $work\$n.log)" }
    } finally { Pop-Location }
    $log = Get-Content -LiteralPath (Join-Path $work "$n.log") -Raw
    if ($log -match '(?i)Missing character|Overfull \\[hv]box|LaTeX Warning|Font shape .* undefined') { throw "Log check failed for $n ($work\$n.log)" }
    Copy-Item -LiteralPath (Join-Path $work "$n.pdf") -Destination (Join-Path $dest "$n.pdf") -Force
    Write-Output "IAM_F02_PREVIEW $n"
}
Write-Output 'IAM_F02_BUILD_DONE'
