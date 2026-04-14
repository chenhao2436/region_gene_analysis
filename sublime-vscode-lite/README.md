# Sublime Text VS Code Lite

这套配置是按你现在的习惯做的轻量版映射：

- 打开项目文件夹后能看到左侧目录树
- `Ctrl+P` 快速切文件
- `Ctrl+B` / `Ctrl+Shift+E` 切换侧边栏
- `Ctrl+N` 新建文件
- `Ctrl+Shift+F` 项目内搜索
- 保留深色风格、等宽字体、自动去尾空格、行尾补换行

## 你要怎么用

1. 在 Sublime Text 里打开 `BSA.sublime-project`
2. 如果你希望更像 VS Code，先把这份配置里的内容复制到 Sublime 的 `User` 配置目录
3. 以后直接打开这个 project 文件，就会按这个工作区启动

## 和 VS Code 的对应关系

- VS Code 的 Explorer -> Sublime 左侧 Side Bar
- VS Code 的 `Ctrl+P` -> Sublime 默认也是 `Ctrl+P`
- VS Code 的 `Ctrl+Shift+P` -> Sublime 命令面板
- VS Code 的多标签切换 -> Sublime 的标签页和 `Ctrl+Tab`

## 我刻意没加的东西

- 没加重型 AI 插件
- 没加 LSP / 补全增强
- 没加复杂构建系统

这样会更轻，也更接近你说的“以查看和简单修改为主”
