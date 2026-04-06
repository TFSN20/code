# 代码示例
## 4个按钮 异步同步 注册/撤销全局快捷键
```
import sys
import time
import ctypes
from PyQt6.QtWidgets import QApplication, QWidget, QVBoxLayout, QPushButton
from PyQt6.QtCore import QThread, pyqtSignal, QAbstractNativeEventFilter

# Win32 常量（底层系统调用参数）
MOD_CONTROL = 0x0002
MOD_ALT = 0x0001
WM_HOTKEY = 0x0312
HOTKEY_ID = 1001

class AsyncWorker(QThread):
    """子线程工作器"""
    finished = pyqtSignal()
    def run(self):
        # 此方法在新 OS 线程中执行
        print("[Async] 开始执行非阻塞任务 (子线程)")
        time.sleep(3)
        print("[Async] 任务结束，发射 finished 信号")
        self.finished.emit()

import ctypes
from ctypes import wintypes
from PyQt6.QtCore import QAbstractNativeEventFilter

# 1. 定义与 Win32 API 完全一致的 MSG 结构体内存布局
class MSG(ctypes.Structure):
    _fields_ = [
        ('hwnd', wintypes.HWND),    # 窗口句柄 (8字节/4字节)
        ('message', wintypes.UINT), # 消息类型 (WM_HOTKEY = 0x0312)
        ('wParam', wintypes.WPARAM),# 附加参数1 (存放热键ID)
        ('lParam', wintypes.LPARAM),# 附加参数2 (存放修饰键与虚拟键码)
        ('time', wintypes.DWORD),   # 消息生成时间戳
        ('pt', wintypes.POINT),     # 鼠标坐标 (未使用)
    ]

class NativeHotkeyFilter(QAbstractNativeEventFilter):
    def __init__(self, callback):
        super().__init__()
        self.callback = callback

    def nativeEventFilter(self, eventType, message):
        if eventType == "windows_generic_MSG":
            # 2. 将 sip.voidptr (裸指针地址) 转换为 ctypes 结构体指针
            msg_ptr = ctypes.cast(int(message), ctypes.POINTER(MSG))
            msg = msg_ptr.contents  # 解引用：读取该内存地址的实际结构体

            # 3. 直接访问底层字段（零字典开销，纯内存读取）
            if msg.message == WM_HOTKEY and msg.wParam == HOTKEY_ID:
                self.callback()
                return True, 0  # (已拦截, Win32 处理返回值)
                
        return False, 0  # 放行给 Qt 默认事件系统

class MainWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("底层原理验证窗口")
        self.resize(420, 300)

        layout = QVBoxLayout(self)
        self.btn_sync = QPushButton("1. 同步打印（阻塞主线程）")
        self.btn_async = QPushButton("2. 非阻塞打印（子线程）")
        self.btn_reg = QPushButton("3. 注册全局快捷键 (Ctrl+Alt+O)")
        self.btn_unreg = QPushButton("4. 撤销全局快捷键")
        
        layout.addWidget(self.btn_sync)
        layout.addWidget(self.btn_async)
        layout.addWidget(self.btn_reg)
        layout.addWidget(self.btn_unreg)

        self.worker = AsyncWorker()
        self.hotkey_filter = None
        self.is_registered = False

        self.btn_sync.clicked.connect(self.on_sync_click)
        self.btn_async.clicked.connect(self.on_async_click)
        self.btn_reg.clicked.connect(self.register_hotkey)
        self.btn_unreg.clicked.connect(self.unregister_hotkey)
        self.worker.finished.connect(lambda: print("[Async] 主线程收到 finished 信号"))

    def on_sync_click(self):
        print("[Sync] 开始执行阻塞任务 (主线程)")
        time.sleep(3)  # 强制挂起当前线程
        print("[Sync] 阻塞结束")

    def on_async_click(self):
        if self.worker.isRunning():
            print("[Async] 任务已在运行")
            return
        self.worker.start()  # 请求 OS 创建新线程并执行 run()

    def register_hotkey(self):
        if self.is_registered:
            return
        # 1. 安装原生事件过滤器
        self.hotkey_filter = NativeHotkeyFilter(self.on_hotkey_triggered)
        QApplication.instance().installNativeEventFilter(self.hotkey_filter)
        
        # 2. 调用 Win32 API 向系统注册全局热键
        # hWnd=None 表示绑定到当前线程消息队列，不依赖窗口焦点
        ret = ctypes.windll.user32.RegisterHotKey(None, HOTKEY_ID, MOD_CONTROL | MOD_ALT, ord('O'))
        if ret:
            self.is_registered = True
            print("[Hotkey] 全局快捷键注册成功 (Ctrl+Alt+O)")
        else:
            print("[Hotkey] 注册失败 (可能热键冲突)")

    def unregister_hotkey(self):
        if not self.is_registered:
            return
        ctypes.windll.user32.UnregisterHotKey(None, HOTKEY_ID)
        QApplication.instance().removeNativeEventFilter(self.hotkey_filter)
        self.hotkey_filter = None
        self.is_registered = False
        print("[Hotkey] 全局快捷键已撤销")

    def on_hotkey_triggered(self):
        print("[Hotkey] 全局快捷键触发！打印成功。")
        self.raise_()  # 触发后将窗口提到前台

if __name__ == "__main__":
    app = QApplication(sys.argv)
    win = MainWindow()
    win.show()
    sys.exit(app.exec())
```
