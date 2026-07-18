# Colab CLI 使用说明(通用,适用于任何项目)

维护:2026-07-19。基于 ROUND63 S2 campaign 实战验证过的全套流程(5 会话
×30 车道零丢失跑完)。所有命令在 **WSL** 里执行;Windows 侧只负责把脚本
文件写到 D 盘再喂给 `wsl.exe`。

---

## 1. 账号与使用策略

| 账号 | 档位 | 并发 GPU 会话 | HOME(CLI 状态目录) | 用途策略 |
|---|---|---|---|---|
| pro1 | Pro | 2 × L4 | `/var/tmp/codex-colab-accounts/pro1` | 正式 campaign 偶数分片 |
| pro2 | **Pro+** | 3 × L4 | `/var/tmp/codex-colab-accounts/pro2` | 正式奇数分片;**自由用途优先烧这个**(算力单元先到期) |

- 正式 campaign 的分片↔账号奇偶对应是冻结约定;随手跑的东西不受限。
- **网页版 Colab 和 CLI 可以并用**(同账号),但不要动
  `/var/tmp/codex-colab-accounts/*` 里的状态文件。

## 2. CLI 基本操作

CLI 本体:`/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab`。
每条命令都要带 `HOME=<账号目录>` 前缀和 `--auth oauth2`,建议再套
`timeout 240`:

```bash
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
H=/var/tmp/codex-colab-accounts/pro2   # 或 pro1

# 新建会话(名字自取,GPU 用 L4)
HOME=$H timeout 240 $C --auth oauth2 new --session my_sess --gpu L4

# 列出会话
HOME=$H timeout 240 $C --auth oauth2 sessions

# 上传 / 下载(本地路径用 /mnt/d/...)
HOME=$H timeout 240 $C --auth oauth2 upload   --session my_sess /mnt/d/work/x.py /content/x.py
HOME=$H timeout 240 $C --auth oauth2 download --session my_sess /content/out.csv /mnt/d/work/out.csv

# 执行(--file 里是发到 VM 上跑的 python;长任务给足 --timeout 秒数)
HOME=$H timeout 300 $C --auth oauth2 exec --session my_sess --file /mnt/d/work/job.py --timeout 14400

# 保活守护(必须 setsid,见 §4)
setsid HOME=$H $C --auth oauth2 keep-alive <ENDPOINT> my_sess >/dev/null 2>&1 &
```

VM 环境:预装 torch/CUDA;**numpy 钉在 1.26.4、scipy 1.13.1**(别在 VM
上乱升级,会撞其它包)。

## 3. 两条铁律(血泪教训)

### 3.1 永远用"文件承载"的 bash 脚本,不要内联单行命令
PowerShell 会把 `wsl.exe -- bash -c "... $VAR ..."` 里的 `$` 变量吃掉
(踩过两次)。正确姿势:

1. 把 bash 脚本写成 D 盘上的 `.sh` 文件;
2. 先 `sed -i 's/\r$//' <script>`(去 CRLF,否则 bash 报错);
3. `wsl.exe -- bash <script>`。

### 3.2 runtime-proxy token 约 30–45 分钟过期,且 CLI 会"自毁"
token 过期后 CLI 收到 401 时**不是刷新而是把会话行从本地状态里剪掉**
(表现:`sessions` 里会话凭空消失、upload/download 404)。VM 本身没死,
计算不受影响——**千万不要因此重建会话**。修复原理(现成脚本
`code/round63/colab/live_rebind.sh`,可直接搬去别的项目):

- 用 `state.client.list_assignments()` 拿当前 assignment(其中
  `runtime_proxy_info` 携带新鲜 token);
- 用 `state.store.add(SessionState(...))` 把会话行写回本地状态;
- 注意 `AssignmentVariant` 是 int 枚举,写回时用 `.name`。

长跑任务的看门狗应当**每个周期先跑一遍 rebind 再查状态**(自愈式,见
`live_watch.sh`)。

## 4. 长跑任务的完整纪律

1. **打包上传**:代码+数据打成 tar 上传解包(大文件分块),见
   `make_bundle.sh`。
2. **驱动脚本在 VM 上自治**:任务清单写在 VM 本地,driver 逐格跑、写
   心跳 JSON + 每格结果落盘(`session_driver.sh` + `remote_lane.py`
   模式),这样 token/网络断了计算也不停。
3. **keep-alive 必须 `setsid`**:WSL 里裸 `nohup` 会随终端死掉;从
   PowerShell 侧启动则用 `Start-Process wsl.exe`。
4. **看门狗 + 定期摘要**:每 15 分钟拉一次心跳 JSON,>10 分钟无心跳报
   停滞;每小时出一行摘要。模板 `live_watch.sh`。
5. **取结果幂等**:按"已完成清单"逐个下载、跳过已存在文件
   (`live_fetch_all.sh` 模式);结果目录结构注意 `shards/<id>.csv` +
   `<id>_meta.json` 命名。
6. **用完释放**:`state.client.unassign(endpoint)`(模板
   `live_release_pro2.sh`),再杀对应 keep-alive。不释放会一直烧算力
   单元。

## 5. 现成脚本地图(`code/round63/colab/`,全部可当模板抄)

| 脚本 | 作用 |
|---|---|
| `make_bundle.sh` | 打包代码+数据为上传 bundle |
| `live_create_sessions.sh` | 批量建 5 会话 |
| `live_launch_all.sh` | 上传 bundle + 启动全部车道 driver |
| `session_driver.sh` / `remote_lane.py` | VM 端自治驱动(心跳+逐格落盘) |
| `live_rebind.sh` | **token 过期自愈**(§3.2) |
| `live_keepalive_attach.sh` | 补挂 keep-alive 守护 |
| `live_watch.sh` | 自愈式看门狗(每周期先 rebind) |
| `live_status_once.sh` | 手动一次性查全会话进度 |
| `live_fetch_all.sh` | 幂等取回全部已完成分片 |
| `live_release_pro2.sh` | 释放 VM + 杀 keep-alive(改 endpoint 列表即可复用) |

## 6. 常见故障速查

| 症状 | 原因 | 处置 |
|---|---|---|
| `sessions` 里会话消失 / 404 | token 过期被 CLI 剪行 | 跑 `live_rebind.sh`,**不要重建会话** |
| upload/download 挂住 | 没套 `timeout`;或 token 过期 | 加 `timeout 240`;先 rebind |
| bash 报 `$'\r'` 错 | CRLF | `sed -i 's/\r$//'` |
| 变量被展开成空 | PowerShell 内联 `$` | 文件承载脚本(§3.1) |
| keep-alive 死了 | WSL nohup 不持久 | `setsid`(§4.3) |
| VM 上 pip 装包冲突 | numpy/scipy 钉版本 | 装包时 `--no-deps` 或换兼容版本 |

## 7. 当前状态快照(2026-07-19 凌晨,过期作废)

- **pro2:全空**(3 槽已释放,assignments: NONE)→ 可随意用。
- **pro1:2 槽被 ROUND63 S2 尾巴占用**,预计柏林时间 03:00–05:00 释放。
