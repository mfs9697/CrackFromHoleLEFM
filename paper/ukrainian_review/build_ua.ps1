# Run from any directory with XeLaTeX, BibTeX, STIX Two Text, STIX Two Math, unicode-math and Ukrainian language support installed.
$ErrorActionPreference = 'Stop'

function Invoke-BuildCommand {
    param([string]$Command, [string[]]$BuildArguments)
    & $Command @BuildArguments
    if ($LASTEXITCODE -ne 0) {
        throw "$Command returned exit code $LASTEXITCODE"
    }
}

Push-Location $PSScriptRoot
try {
    $job = 'AAM_UA_Bogdanov_review'
    $work = Join-Path $PSScriptRoot '.build'
    $null = New-Item -ItemType Directory -Path $work -Force
    $xeArgs = @('-interaction=nonstopmode', '-halt-on-error', '-file-line-error',
                '-recorder', "-output-directory=$work", "-jobname=$job", 'main_ua.tex')
    Invoke-BuildCommand xelatex $xeArgs

    $previousBibInputs = $env:BIBINPUTS
    Push-Location $work
    try {
        $env:BIBINPUTS = "$PSScriptRoot$([IO.Path]::PathSeparator)$previousBibInputs"
        Invoke-BuildCommand bibtex @($job)
    }
    finally {
        $env:BIBINPUTS = $previousBibInputs
        Pop-Location
    }

    Invoke-BuildCommand xelatex $xeArgs
    Invoke-BuildCommand xelatex $xeArgs
    $log = Get-Content -LiteralPath (Join-Path $work "$job.log") -Raw
    if ($log -match '(?i)undefined (references|citations)|multiply defined|Overfull \\[hv]box|Missing character|No hyphenation patterns') {
        throw "The final LaTeX log contains unresolved references, overflow, or a font/language issue: $work"
    }
    Copy-Item -LiteralPath (Join-Path $work "$job.pdf") -Destination (Join-Path $PSScriptRoot "$job.pdf") -Force
}
finally {
    Pop-Location
}
