
#pragma once

////                Window               ////

void SetWindowSize(int width, int height);
void SetWindowPosition(int x, int y);
void SetWindowResizeCallback(void(*callback)(int, int));
void SetWindowMoveCallback(void(*callback)(int, int));
void SetWindowName(const char* name);

void GetWindowSize(int* x, int* y);
void GetWindowPos(int* x, int* y);

// Sets Window as Full Screen, with given resolution, this might improve performance but pixel density drops down
// You can set the resolution using GetMonitorSize or following resolutions: 
// 2560x1440, 1920x1080, 1280×720 etc.
bool EnterFullscreen(int fullscreenWidth, int fullscreenHeight);

// Go back to full screen Mode
bool ExitFullscreen(int windowX, int windowY, int windowedWidth, int windowedHeight);

void GetMonitorSize(int* width, int* height);
void SetFocusChangedCallback(void(*callback)(bool));

void TerminateWindow();
void HandleInput();
////                Keyboard             ////

// use KeyboardKey enum or asci char 'X'
bool GetKeyDown(char c);
// use KeyboardKey enum or asci char 'X'
bool GetKeyPressed(char c);
// use KeyboardKey enum or asci char 'X'
bool GetKeyReleased(char c);

void SetKeyPressCallback(void(*callback)(wchar_t));

////                Mouse                ////

enum MouseButton_
{
    MouseButton_Left   = 1,
    MouseButton_Right  = 2,
    MouseButton_Middle = 4
};
typedef int MouseButton;

bool GetMouseDown(MouseButton button);
bool GetMouseReleased(MouseButton button);
bool GetMousePressed(MouseButton button);

void SetMousePos(float x, float y);
void GetMousePos(float* x, float* y);
void GetMouseWindowPos(float* x, float* y);
void SetMouseWindowPos(float x, float y);
void SetMouseMoveCallback(void(*callback)(float, float));

////                TIME                 ////

float GetDeltaTime();

extern struct android_app* g_android_app;

/*
* VK_0 - VK_9 are the same as ASCII '0' - '9' (0x30 - 0x39)
* 0x3A - 0x40 : unassigned
* VK_A - VK_Z are the same as ASCII 'A' - 'Z' (0x41 - 0x5A)
*/

// Todo: we will need different enum values for each platform

enum KeyboardKey_
{
    Key_BACK           = 0x08,
    Key_TAB            = 0x09,
    Key_CLEAR          = 0x0C,
    Key_RETURN         = 0x0D,
    Key_SHIFT          = 0x10,
    Key_CONTROL        = 0x11,
    Key_MENU           = 0x12,
    Key_PAUSE          = 0x13,
    Key_CAPITAL        = 0x14,
    Key_ESCAPE         = 0x1B,
    Key_CONVERT        = 0x1C,
    Key_NONCONVERT     = 0x1D,
    Key_ACCEPT         = 0x1E,
    Key_MODECHANGE     = 0x1F,
    Key_SPACE          = 0x20,
    Key_PRIOR          = 0x21,
    Key_NEXT           = 0x22,
    Key_END            = 0x23,
    Key_HOME           = 0x24,
    Key_LEFT           = 0x25,
    Key_UP             = 0x26,
    Key_RIGHT          = 0x27,
    Key_DOWN           = 0x28,
    Key_SELECT         = 0x29,
    Key_PRINT          = 0x2A,
    Key_EXECUTE        = 0x2B,
    Key_SNAPSHOT       = 0x2C,
    Key_INSERT         = 0x2D,
    Key_DELETE         = 0x2E,
    Key_HELP           = 0x2F,
    Key_LWIN           = 0x5B,
    Key_RWIN           = 0x5C,
    Key_APPS           = 0x5D,
    Key_SLEEP          = 0x5F,
    Key_NUMPAD0        = 0x60,
    Key_NUMPAD1        = 0x61,
    Key_NUMPAD2        = 0x62,
    Key_NUMPAD3        = 0x63,
    Key_NUMPAD4        = 0x64,
    Key_NUMPAD5        = 0x65,
    Key_NUMPAD6        = 0x66,
    Key_NUMPAD7        = 0x67,
    Key_NUMPAD8        = 0x68,
    Key_NUMPAD9        = 0x69,
    Key_MULTIPLY       = 0x6A,
    Key_ADD            = 0x6B,
    Key_SEPARATOR      = 0x6C,
    Key_SUBTRACT       = 0x6D,
    Key_DECIMAL        = 0x6E,
    Key_DIVIDE         = 0x6F,
    Key_F1             = 0x70,
    Key_F2             = 0x71,
    Key_F3             = 0x72,
    Key_F4             = 0x73,
    Key_F5             = 0x74,
    Key_F6             = 0x75,
    Key_F7             = 0x76,
    Key_F8             = 0x77,
    Key_F9             = 0x78,
    Key_F10            = 0x79,
    Key_F11            = 0x7A,
    Key_F12            = 0x7B,
};
typedef int KeyboardKey;