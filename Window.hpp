
#pragma once

void SetWindowSize(int width, int height);

void SetWindowPosition(int x, int y);

void SetWindowResizeCallback(void(*callback)(int, int));

void SetWindowMoveCallback(void(*callback)(int, int));

void GetWindowSize(int* x, int* y);

void GetWindowPos(int* x, int* y);

void SetWindowName(const char* name);

void TerminateWindow();

extern struct android_app* g_android_app;