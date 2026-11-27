# Build the Ukrainian figure panels (AAM-UA7) from figures_src/panels/**/*.tex with XeLaTeX.
# Each panel is compiled in its own folder into figures_src/.build/<group>/ and copied to figures/<group>/<name>.pdf.
# Requires XeLaTeX with pgfplots, fontspec, unicode-math and the STIX Two OpenType fonts (stix2-otf).
# Usage (repository root or any directory):
#   .\paper\ukrainian_review\scripts\build_ua_figures.ps1                  # all panels
#   .\paper\ukrainian_review\scripts\build_ua_figures.ps1 -Only trajectory/late_detail
#   .\paper\ukrainian_review\scripts\build_ua_figures.ps1 -SourceRoot <copy> -OutputRoot <dir>   (negative controls)
param([string[]]$Only = @(), [string]$SourceRoot = '', [string]$OutputRoot = '')
$ErrorActionPreference = 'Stop'
$ua = Split-Path -Parent $PSScriptRoot
if (-not $SourceRoot) { $SourceRoot = Join-Path $ua 'figures_src' }
if (-not $OutputRoot) { $OutputRoot = Join-Path $ua 'figures' }
$panels = Get-ChildItem -LiteralPath (Join-Path $SourceRoot 'panels') -Recurse -Filter *.tex | Sort-Object FullName
$built = 0
foreach ($p in $panels) {
    $group = Split-Path -Leaf $p.DirectoryName
    $id = "$group/$($p.BaseName)"
    if ($Only.Count -gt 0 -and -not ($Only -contains $id)) { continue }
    $work = Join-Path $SourceRoot (Join-Path '.build' $group)
    $null = New-Item -ItemType Directory -Path $work -Force
    Push-Location $p.DirectoryName
    try {
        & xelatex -interaction=nonstopmode -halt-on-error -file-line-error "-output-directory=$work" $p.Name | Out-Null
        if ($LASTEXITCODE -ne 0) { throw "xelatex failed for $id (see $work\$($p.BaseName).log)" }
    }
    finally { Pop-Location }
    $log = Get-Content -LiteralPath (Join-Path $work "$($p.BaseName).log") -Raw
    if ($log -match '(?i)Missing character|Overfull \\hbox|Package pgfplots Warning|Font shape .* undefined|LaTeX Warning') {
        throw "Log check failed for ${id}: missing glyph, overflow, pgfplots or LaTeX warning ($work\$($p.BaseName).log)"
    }
    $dest = Join-Path $OutputRoot $group
    $null = New-Item -ItemType Directory -Path $dest -Force
    Copy-Item -LiteralPath (Join-Path $work "$($p.BaseName).pdf") -Destination (Join-Path $dest "$($p.BaseName).pdf") -Force
    $built++
    Write-Output "UA7_PANEL $id"
}
Write-Output "UA7_PANELS_BUILT $built"
