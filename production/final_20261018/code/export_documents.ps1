param([string]$RunRoot = (Split-Path -Parent $PSScriptRoot))
$ErrorActionPreference = 'Stop'
$env:SystemRoot = 'C:\Windows'
$env:windir = 'C:\Windows'
$env:ComSpec = 'C:\Windows\System32\cmd.exe'
$env:LOCALAPPDATA = Join-Path $env:USERPROFILE 'AppData\Local'
$env:APPDATA = Join-Path $env:USERPROFILE 'AppData\Roaming'
$env:TEMP = Join-Path $env:LOCALAPPDATA 'Temp'
$out = if ($env:ICI_SUBMISSION_ROOT) { $env:ICI_SUBMISSION_ROOT } else { Join-Path $RunRoot 'submission_final' }
$dest = Join-Path $RunRoot 'qa\word'
New-Item -ItemType Directory -Force -Path $dest | Out-Null
$word = New-Object -ComObject Word.Application
$word.Visible = $false
$word.DisplayAlerts = 0
$records = @()
try {
 foreach ($name in @('Manuscript','Cover_Letter','Table_1')) {
  $inputPath = Join-Path $out ($name+'.docx')
  $outputPath = Join-Path $dest ($name+'.pdf')
  $started = (Get-Date).ToUniversalTime().ToString('o')
  $doc = $word.Documents.Open($inputPath,$false,$true)
  try {
   $doc.Repaginate()
   $doc.ExportAsFixedFormat($outputPath,17)
   $records += [pscustomobject]@{input=$inputPath;output=$outputPath;renderer='Microsoft Word';version=$word.Version;started=$started;finished=(Get-Date).ToUniversalTime().ToString('o');status='EXPORTED'}
  } finally { $doc.Close(0) }
 }
} finally { $word.Quit() }
$records | ConvertTo-Json -Depth 4 | Set-Content -LiteralPath (Join-Path $RunRoot 'qa\word_export_record.json') -Encoding UTF8
Write-Output 'Word PDF exports complete'
