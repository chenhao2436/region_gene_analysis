# Open In Sublime

Use this launcher when you want a generic way to open any current project or file in the installed Sublime Text.

## Usage

Run it with no arguments to open the current working directory:

```bat
open-in-sublime.cmd
```

Run it with a folder or file path to open that target directly:

```bat
open-in-sublime.cmd D:\app\.codex\project\BSA
open-in-sublime.cmd D:\app\.codex\project\BSA\ploidy30_bsa_bundle_20260318\README.md
```

## Notes

- This always targets the installed Sublime Text at `D:\app\Sublime Text\sublime_text.exe`.
- If you run it from a Codex task or run window without arguments, it opens the current directory.
