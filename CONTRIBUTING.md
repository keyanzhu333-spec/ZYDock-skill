# Contributing to ZYDock / 贡献指南

Thanks for your interest in improving ZYDock! Contributions of all kinds are welcome — bug reports, new scripts, documentation fixes, or ideas.

感谢你有兴趣改进 ZYDock!欢迎各种形式的贡献 —— 报告 bug、新增脚本、修正文档或提出想法。

---

## Reporting bugs & requesting features / 报告问题与提需求

Open an [Issue](https://github.com/keyanzhu333-spec/ZYDock-skill/issues) and include:

请在 [Issues](https://github.com/keyanzhu333-spec/ZYDock-skill/issues) 里提交,并尽量包含:

- What you tried (the command you ran) / 你执行的命令
- What you expected vs. what happened / 期望的结果 vs. 实际的结果
- Error messages / logs (redact any API keys!) / 报错信息或日志(**记得抹掉 API key!**)
- Your environment: OS, Python version / 运行环境:操作系统、Python 版本

---

## Submitting code (Pull Request) / 提交代码(Pull Request)

You do **not** need write access. Use the standard fork + PR flow:

你**不需要**仓库写权限,按标准的 fork + PR 流程即可:

```bash
# 1. Fork this repo on GitHub (top-right "Fork" button)
#    在 GitHub 上点右上角 "Fork" 按钮

# 2. Clone your fork / 克隆你自己的 fork
git clone https://github.com/<your-username>/ZYDock-skill.git
cd ZYDock-skill

# 3. Create a branch — don't work on main / 新建分支,不要直接改 main
git checkout -b fix-something

# 4. Make changes, then commit / 修改后提交
git add .
git commit -m "Describe what you changed"

# 5. Push to your fork / 推送到你的 fork
git push origin fix-something

# 6. Open a Pull Request on GitHub / 在 GitHub 上发起 Pull Request
```

The maintainer will review the diff, may leave comments, and merge once it looks good.

维护者会审阅改动,可能会留言讨论,确认没问题后合并。

---

## Guidelines / 约定

- **One PR = one topic.** Keep unrelated changes in separate PRs. / 一个 PR 只做一件事,无关改动请拆开。
- **Never commit secrets.** No API keys, credentials, or `.env` files — they're already in `.gitignore`. Real keys live outside the repo in `~/.config/ZYDock/`. / **绝不提交密钥。** 不要提交 API key、凭据或 `.env`(已在 `.gitignore` 中);真实密钥存放在仓库之外的 `~/.config/ZYDock/`。
- **Match the existing style.** Scripts are Python 3, use `argparse` for CLI args, and follow the patterns in [`scripts/`](scripts/). / 保持风格一致:脚本用 Python 3,命令行参数用 `argparse`,参考 [`scripts/`](scripts/) 里的现有写法。
- **Test before submitting.** Run your changed script end-to-end and confirm it works. / 提交前请端到端跑一遍你改动的脚本,确认可用。
- **Update docs.** If you add or change a script, update the relevant table/section in [`README.md`](README.md) and [`SKILL.md`](SKILL.md). / 若新增或修改脚本,请同步更新 [`README.md`](README.md) 和 [`SKILL.md`](SKILL.md) 里的相关表格/说明。

---

## License / 许可协议

By contributing, you agree that your contributions will be licensed under the [MIT License](LICENSE) that covers this project.

提交贡献即表示你同意你的贡献以本项目采用的 [MIT 许可协议](LICENSE)授权。
