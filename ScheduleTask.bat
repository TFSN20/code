@echo off
:: ����Ƿ��Թ���Ա�������
net session >nul 2>&1
if %errorlevel% neq 0 (
    echo ��ȴ�...
    powershell start-process -filepath '%0' -verb runas
    exit /b
)

set taskname=ConnectWiFi
set taskpath=C:\Users\Administrator\Desktop\connect-wifi.bat
set time=06:59

schtasks /create /tn "%taskname%" /tr "%taskpath%" /sc daily /st %time% /f

set message=���� "%taskname%" �ѳɹ�����������ÿ������%time%���С�
echo %message%
pause
