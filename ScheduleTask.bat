@echo off
:: 检查是否以管理员身份运行
net session >nul 2>&1
if %errorlevel% neq 0 (
    echo 请等待...
    powershell start-process -filepath '%0' -verb runas
    exit /b
)

set taskname=ConnectWiFi
set taskpath=C:\Users\Administrator\Desktop\connect-wifi.bat
set time=06:59

schtasks /create /tn "%taskname%" /tr "%taskpath%" /sc daily /st %time% /f

set message=任务 "%taskname%" 已成功创建，并将每天上午%time%运行。
echo %message%
pause
