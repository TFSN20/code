## 版本
- Auto.js_4.1.1 Alpha2：支持importClass，不支持模板字符串
- Autojsx：最新版不支持importClass，支持模板字符串。
## 常见问题
- http模块请求时最好不携带accept-encoding请求头，否则响应体大概率乱码。
- 关于悬浮窗和UI上的界面操作（按钮点击），如果有耗时操作会堵塞，使其不能点击，使用threads.start(()=>{})解决。
- jsBridge.callHandler里的网络请求必须threads.start。
- 悬浮窗按钮点击后触发的函数若有堵塞，则必须先 yield myAwait(sleep_(10));
