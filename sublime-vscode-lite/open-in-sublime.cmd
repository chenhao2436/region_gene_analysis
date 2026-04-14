@echo off
setlocal

set "SUBLIME_EXE=D:\app\Sublime Text\sublime_text.exe"
set "TARGET=%~1"

if not exist "%SUBLIME_EXE%" (
  echo Sublime Text not found: %SUBLIME_EXE%
  exit /b 1
)

if "%TARGET%"=="" (
  start "" "%SUBLIME_EXE%" "%CD%"
) else (
  if exist "%TARGET%\NUL" (
    start "" "%SUBLIME_EXE%" "%TARGET%"
  ) else (
    start "" "%SUBLIME_EXE%" "%TARGET%"
  )
)
