
#ifdef __ANDROID__
    #include <game-activity/native_app_glue/android_native_app_glue.h>
    #include <EGL/egl.h>
    #include <GLES3/gl32.h>
    #define STBI_NO_STDIO
    #define STBI_NEON

    static EGLDisplay display_ = EGL_NO_DISPLAY;
    static EGLSurface surface_ = EGL_NO_SURFACE;
    static EGLContext context_ = EGL_NO_CONTEXT;
    android_app* g_android_app = nullptr;
#else
    #include "glad.hpp"

    #define AX_LOG(...) printf(__VA_ARGS__)
#endif

#include "Memory.hpp"
#include "Renderer.hpp"
#include "Window.hpp"
#include <stdlib.h> // exit failure

int windowWidth_   = 1240;
int windowHeight_  = 720;
int windowPosX_    = 500;
int windowPosY_    = 500;

static char WindowName[128]{ 'A', 'S', 'T', 'L' };
static void(*WindowMoveCallback)  (int, int) = nullptr;
static void(*WindowResizeCallback)(int, int) = nullptr;
static void(*MouseMoveCallback)(float, float) = nullptr;
static void(*KeyPressCallback)(wchar_t) = nullptr;
static void(*FocusChangedCallback)(bool) = nullptr;

extern void AXInit();
extern int  AXStart();
extern void AXLoop();
extern void AXExit();

void UpdateRenderArea();

#define AX_USE_WINDOW

#if defined(AX_USE_WINDOW) && defined(__ANDROID__) 
static void InitWindow()
{
  constexpr EGLint attribs[] = {
                               EGL_RENDERABLE_TYPE, EGL_OPENGL_ES3_BIT,
                               EGL_SURFACE_TYPE, EGL_WINDOW_BIT,
                               EGL_BLUE_SIZE, 8,
                               EGL_GREEN_SIZE, 8,
                               EGL_RED_SIZE, 8,
                               EGL_DEPTH_SIZE, 24,
                               EGL_NONE
  };

  // The default display is probably what you want on Android
  EGLDisplay display = eglGetDisplay(EGL_DEFAULT_DISPLAY);
  eglInitialize(display, nullptr, nullptr);

  // figure out how many configs there are
  EGLint numConfigs;
  eglChooseConfig(display, attribs, nullptr, 0, &numConfigs);

  // get the list of configurations
  EGLConfig supportedConfigs[32]{};
  eglChooseConfig(display, attribs, supportedConfigs, numConfigs, &numConfigs);

  // Find a config we like.
  // Could likely just grab the first if we don't care about anything else in the config.
  // Otherwise hook in your own heuristic
  EGLConfig config = nullptr;
  for (int i = 0; i < numConfigs; i++)
  {
    EGLint red, green, blue, depth;

    config = supportedConfigs[i];
    if (eglGetConfigAttrib(display, config, EGL_RED_SIZE, &red)
      && eglGetConfigAttrib(display, config, EGL_GREEN_SIZE, &green)
      && eglGetConfigAttrib(display, config, EGL_BLUE_SIZE, &blue)
      && eglGetConfigAttrib(display, config, EGL_DEPTH_SIZE, &depth))
      if (red == 8 && green == 8 && blue == 8 && depth == 24)
        break;
  }
  AX_LOG("Found %i configs\n", numConfigs);

  // create the proper window surface
  EGLint format;
  eglGetConfigAttrib(display, config, EGL_NATIVE_VISUAL_ID, &format);
  EGLSurface surface = eglCreateWindowSurface(display, config, g_android_app->window, nullptr);

  // Create a GLES 3 context
  EGLint contextAttribs[] = { EGL_CONTEXT_CLIENT_VERSION, 3, EGL_NONE };
  EGLContext context = eglCreateContext(display, config, nullptr, contextAttribs);

  // get some window metrics
  EGLBoolean madeCurrent = eglMakeCurrent(display, surface, surface, context);

  display_ = display;
  surface_ = surface;
  context_ = context;

  // make width and height invalid so it gets updated the first frame in @a updateRenderArea()
  windowWidth_ = -1;
  windowHeight_ = -1;
}

void SetWindowSize(int width, int height) { }

void SetWindowPosition(int x, int y) { }

void SetWindowResizeCallback(void(*callback)(int, int)) {}

void SetWindowMoveCallback(void(*callback)(int, int)) {}

void SetMouseMoveCallback(void(*callback)(float, float)) {}

void GetWindowSize(int* x, int* y) { }

void GetWindowPos(int* x, int* y) { }

void SetWindowName(const char* name) { }

#elif defined(AX_USE_WINDOW) && defined(_WIN32)

#define NOMINMAX
#define WIN32_LEAN_AND_MEAN 
#define VC_EXTRALEAN
#include <Windows.h>

#pragma comment (lib, "gdi32.lib")
#pragma comment (lib, "user32.lib")
#pragma comment (lib, "opengl32.lib")

static HWND hwnd = nullptr;

void SetWindowSize(int width, int height)
{
    windowWidth_ = width; windowHeight_ = height;
    if (!hwnd) return;
    SetWindowPos(hwnd, nullptr, windowPosX_, windowPosY_, width, height, 0);
}

void SetWindowPosition(int x, int y)
{
    windowPosX_ = x; windowPosY_ = y;
    if (!hwnd) return;
    SetWindowPos(hwnd, nullptr, x, y, windowWidth_, windowHeight_, 0);
}

void SetFocusChangedCallback(void(*callback)(bool))      { FocusChangedCallback = callback; }
void SetKeyPressCallback(void(*callback)(wchar_t))       { KeyPressCallback     = callback; }
void SetMouseMoveCallback(void(*callback)(float, float)) { MouseMoveCallback    = callback; }
void SetWindowResizeCallback(void(*callback)(int, int))  { WindowResizeCallback = callback; }
void SetWindowMoveCallback(void(*callback)(int, int))    { WindowMoveCallback   = callback; }

void GetWindowSize(int* x, int* y) { *x = windowWidth_; *y = windowHeight_; }
void GetWindowPos(int* x, int* y)  { *x = windowPosX_; *y = windowPosY_;    }

void SetWindowName(const char* name)
{
    SmallMemSet(WindowName, 0, sizeof(WindowName));
    
    for (int i = 0; *name; i++)
    {
        WindowName[i] = *name++;
    }
    if (hwnd)
    {
        SetWindowText(hwnd, WindowName);
    }
}

void GetMonitorSize(int* width, int* height)
{
    *width  = GetSystemMetrics(SM_CXSCREEN);
    *height = GetSystemMetrics(SM_CYSCREEN);
}

// https://stackoverflow.com/questions/2382464/win32-full-screen-and-hiding-taskbar
bool EnterFullscreen(int fullscreenWidth, int fullscreenHeight) 
{
    DEVMODE fullscreenSettings;
    EnumDisplaySettings(NULL, 0, &fullscreenSettings);
    windowWidth_  = fullscreenSettings.dmPelsWidth  = fullscreenWidth;
    windowHeight_ = fullscreenSettings.dmPelsHeight = fullscreenHeight;
    fullscreenSettings.dmFields           = DM_PELSWIDTH | DM_PELSHEIGHT;

    SetWindowLongPtr(hwnd, GWL_EXSTYLE, WS_EX_APPWINDOW | WS_EX_TOPMOST);
    SetWindowLongPtr(hwnd, GWL_STYLE, WS_POPUP | WS_VISIBLE);
    SetWindowPos(hwnd, HWND_TOPMOST, 0, 0, fullscreenWidth, fullscreenHeight, SWP_SHOWWINDOW);
    bool success = ChangeDisplaySettings(&fullscreenSettings, CDS_FULLSCREEN) != DISP_CHANGE_SUCCESSFUL;
    ASSERT(success && "unable to make full screen");
    ShowWindow(hwnd, SW_MAXIMIZE);
    if (success && WindowResizeCallback) WindowResizeCallback(fullscreenWidth, fullscreenHeight);
    return success;
}

bool ExitFullscreen(int windowX, int windowY, int windowedWidth, int windowedHeight) 
{
    SetWindowLongPtr(hwnd, GWL_EXSTYLE, WS_EX_LEFT);
    SetWindowLongPtr(hwnd, GWL_STYLE, WS_OVERLAPPEDWINDOW | WS_VISIBLE);
    bool success = ChangeDisplaySettings(NULL, CDS_RESET) == DISP_CHANGE_SUCCESSFUL;
    ASSERT(success && "unable to make windowed");
    windowWidth_ = windowedWidth; windowHeight_ = windowedHeight;
    SetWindowPos(hwnd, HWND_NOTOPMOST, windowX, windowY, windowedWidth, windowedHeight, SWP_SHOWWINDOW);
    ShowWindow(hwnd, SW_RESTORE);
    if (success && WindowResizeCallback) WindowResizeCallback(windowedWidth, windowedHeight);
    return success;
}

// See https://www.khronos.org/registry/OpenGL/extensions/ARB/WGL_ARB_create_context.txt for all values
// See https://www.khronos.org/registry/OpenGL/extensions/ARB/WGL_ARB_pixel_format.txt for all values
// See https://gist.github.com/nickrolfe/1127313ed1dbf80254b614a721b3ee9c
typedef HGLRC WINAPI wglCreateContextAttribsARB_type(HDC hdc, HGLRC hShareContext, const int* attribList);
wglCreateContextAttribsARB_type* wglCreateContextAttribsARB;

typedef BOOL WINAPI wglChoosePixelFormatARB_type(HDC hdc, const int* piAttribIList, const FLOAT* pfAttribFList, UINT nMaxFormats, int* piFormats, UINT* nNumFormats);
wglChoosePixelFormatARB_type* wglChoosePixelFormatARB = nullptr;

BOOL(WINAPI* wglSwapIntervalEXT)(int) = nullptr;

#include <process.h> // exit
static void FatalError(const char* msg)
{
    MessageBoxA(NULL, msg, "Error", MB_OK | MB_ICONEXCLAMATION);
    exit(EXIT_FAILURE);
}

static void InitOpenGLExtensions(void)
{
    // Before we can load extensions, we need a dummy OpenGL context, created using a dummy window.
    // We use a dummy window because you can only set the pixel format for a window once.
    WNDCLASSA window_class{};
    window_class.style         = CS_HREDRAW | CS_VREDRAW | CS_OWNDC;
    window_class.lpfnWndProc   = DefWindowProcA;
    window_class.hInstance     = GetModuleHandle(0);
    window_class.lpszClassName = "Dummy_WGL_StagingWindow";
    
    if (!RegisterClass(&window_class)) 
        FatalError("Failed to register dummy OpenGL window.");
    
    HWND dummy_window = CreateWindowExA(0, window_class.lpszClassName, "ASTL Window",
                                        0, CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT, 0,
                                        0, window_class.hInstance, 0);
    
    if (!dummy_window)
        FatalError("Failed to create dummy OpenGL window.");
    
    HDC dummy_dc = GetDC(dummy_window);
    PIXELFORMATDESCRIPTOR pfd{};
    pfd.nSize        = sizeof(pfd);
    pfd.nVersion     = 1;
    pfd.iPixelType   = PFD_TYPE_RGBA;
    pfd.dwFlags      = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER;
    pfd.cColorBits   = 32;
    pfd.cAlphaBits   = 8;
    pfd.iLayerType   = PFD_MAIN_PLANE;
    pfd.cDepthBits   = 24;
    pfd.cStencilBits = 8;
    
    int pixel_format = ChoosePixelFormat(dummy_dc, &pfd);
    if (!pixel_format) FatalError("Failed to find a suitable pixel format.");
    
    if (!SetPixelFormat(dummy_dc, pixel_format, &pfd)) FatalError("Failed to set the pixel format.");
    
    HGLRC dummy_context = wglCreateContext(dummy_dc);
    if (!dummy_context) FatalError("Failed to create a dummy OpenGL rendering context.");
    
    if (!wglMakeCurrent(dummy_dc, dummy_context)) FatalError("Failed to activate dummy OpenGL rendering context.");
    
    wglCreateContextAttribsARB = (wglCreateContextAttribsARB_type*)wglGetProcAddress("wglCreateContextAttribsARB");
    wglChoosePixelFormatARB    = (wglChoosePixelFormatARB_type*)wglGetProcAddress("wglChoosePixelFormatARB");
    wglSwapIntervalEXT         = (BOOL(WINAPI*)(int)) wglGetProcAddress("wglSwapIntervalEXT");

    wglMakeCurrent(dummy_dc, 0);
    wglDeleteContext(dummy_context);
    ReleaseDC(dummy_window, dummy_dc);
    DestroyWindow(dummy_window);
}

static HGLRC InitOpenGL(HDC real_dc)
{
    InitOpenGLExtensions();
    // Now we can choose a pixel format the modern way, using wglChoosePixelFormatARB.
    int pixel_format_attribs[] = {
        0x2001,          1, // WGL_DRAW_TO_WINDOW_ARB
        0x2010,          1, // WGL_SUPPORT_OPENGL_ARB
        0x2011,          1, // WGL_DOUBLE_BUFFER_ARB
        0x2003,     0x2027,  // WGL_ACCELERATION_ARB, WGL_FULL_ACCELERATION_ARB
        0x2013,     0x202B,  // WGL_PIXEL_TYPE_ARB, WGL_TYPE_RGBA_ARB
        0x2014,         32,  // WGL_COLOR_BITS_ARB
        0x2022,         24,  // WGL_DEPTH_BITS_ARB
        0x2023,          8,  // WGL_STENCIL_BITS_ARB
        0x20A9,          1,  // WGL_FRAMEBUFFER_SRGB_CAPABLE_ARB <- SRGB support
        0x2041,          1,  // WGL_SAMPLE_BUFFERS_ARB           <- enable MSAA
        0x2042,          8,  // WGL_SAMPLES_ARB                  <- 8x MSAA
        0
    };
    
    int pixel_format;
    UINT num_formats;
    wglChoosePixelFormatARB(real_dc, pixel_format_attribs, 0, 1, &pixel_format, &num_formats);
    if (!num_formats) 
      FatalError("Failed to set the OpenGL 3.3 pixel format.");
    
    PIXELFORMATDESCRIPTOR pfd;
    DescribePixelFormat(real_dc, pixel_format, sizeof(pfd), &pfd);
    if (!SetPixelFormat(real_dc, pixel_format, &pfd)) 
      FatalError("Failed to set the OpenGL 3.3 pixel format.");
    
    // Specify that we want to create an OpenGL 3.2 core profile context
    int gl32_attribs[] = {
        0x2091, 3, // WGL_CONTEXT_MAJOR_VERSION_ARB
        0x2092, 2, // WGL_CONTEXT_MINOR_VERSION_ARB
        0x9126,  0x00000001, // WGL_CONTEXT_PROFILE_MASK_ARB, WGL_CONTEXT_CORE_PROFILE_BIT_ARB
        0,
    };
    
    HGLRC gl32_context = wglCreateContextAttribsARB(real_dc, 0, gl32_attribs);
    if (!gl32_context) 
      FatalError("Failed to create OpenGL 3.2 context.");
    
    if (!wglMakeCurrent(real_dc, gl32_context)) 
      FatalError("Failed to activate OpenGL 3.2 rendering context.");
    
    return gl32_context;
}

// Input Code
unsigned long g_axDownKeys[2]{};
unsigned long g_axLastKeys[2]{};
unsigned long g_axPressedKeys[2]{};
unsigned long g_axReleasedKeys[2]{};
// mouse
int g_axMouseDown = 0, g_axMouseLast = 0, g_axMousePressed = 0, g_axMouseReleased = 0;

inline bool GetBit128(unsigned long bits[2], int idx)   { return !!(bits[idx > 63] & (1ul << (idx & 63ul))); }
inline void SetBit128(unsigned long bits[2], int idx)   { bits[idx > 63] |= 1ul << (idx & 63); }
inline void ResetBit128(unsigned long bits[2], int idx) { bits[idx > 63] &= ~(1ul << (idx & 63)); }

bool GetKeyDown(char c)     { return GetBit128(g_axDownKeys, c); }

bool GetKeyReleased(char c) { return GetBit128(g_axReleasedKeys, c); }

bool GetKeyPressed(char c)  { return GetBit128(g_axPressedKeys, c); }

static void SetPressedAndReleasedKeys()
{
    g_axReleasedKeys[0] = g_axLastKeys[0] & ~g_axDownKeys[0];
    g_axReleasedKeys[1] = g_axLastKeys[1] & ~g_axDownKeys[1];
    g_axPressedKeys[0] = ~g_axLastKeys[0] & g_axDownKeys[0];
    g_axPressedKeys[1] = ~g_axLastKeys[1] & g_axDownKeys[1];
    // Mouse
    g_axMouseReleased = g_axMouseLast & ~g_axMouseDown;
    g_axMousePressed  = ~g_axMouseLast & g_axMouseDown;
}

static void RecordLastKeys() {
    g_axLastKeys[0] = g_axDownKeys[0];
    g_axLastKeys[1] = g_axDownKeys[1];
 
    g_axMouseLast = g_axMouseDown;
}

bool GetMouseDown(MouseButton button)     { return !!(g_axMouseDown & button); }
bool GetMouseReleased(MouseButton button) { return !!(g_axMouseReleased & button); }
bool GetMousePressed(MouseButton button)  { return !!(g_axMousePressed & button); }

void GetMousePos(float* x, float* y)
{
    ASSERT((uint64_t)x & (uint64_t)y); // shouldn't be nullptr
    POINT point;
    GetCursorPos(&point);
    *x = (float)point.x;
    *y = (float)point.y;
}

void SetMousePos(float x, float y)
{
    SetCursorPos((int)x, (int)y);
}

float g_axMousePosX = 0.0, g_axMousePosY = 0.0f;

void GetMouseWindowPos(float* x, float* y)
{
    *x = g_axMousePosX; *y = g_axMousePosY;
}

void SetMouseWindowPos(float x, float y)
{
    SetMousePos(windowPosX_ + x, windowPosY_ + y);
}

static LRESULT CALLBACK WindowCallback(HWND window, UINT msg, WPARAM wparam, LPARAM lparam)
{
    LRESULT result = 0;
    wchar_t wch = 0;

    switch (msg)
    {
        case WM_MOUSEMOVE:
            g_axMousePosX = (float)LOWORD(lparam); 
            g_axMousePosY = (float)HIWORD(lparam); 
            if (MouseMoveCallback) MouseMoveCallback(g_axMousePosX, g_axMousePosY);
            break;
        case WM_LBUTTONDOWN: g_axMouseDown |= MouseButton_Left; break;
        case WM_RBUTTONDOWN: g_axMouseDown |= MouseButton_Right; break;
        case WM_MBUTTONDOWN: g_axMouseDown |= MouseButton_Middle; break;
        case WM_LBUTTONUP:   g_axMouseDown &= ~MouseButton_Left; break;
        case WM_RBUTTONUP:   g_axMouseDown &= ~MouseButton_Right; break;
        case WM_MBUTTONUP:   g_axMouseDown &= ~MouseButton_Middle; break;
        case WM_XBUTTONUP:   g_axMouseDown &= ~MouseButton_Middle; break;
        case WM_KEYDOWN:
        case WM_SYSKEYDOWN:
        {
            if (wparam > 127) break;
            SetBit128(g_axDownKeys, wparam);
            break;
        }
        case WM_SETFOCUS:
        case WM_KILLFOCUS:
            if (FocusChangedCallback) FocusChangedCallback(msg == WM_SETFOCUS);
            break;
        case WM_KEYUP:
        case WM_SYSKEYUP:
        {
            if (wparam > 127) break;
            ResetBit128(g_axDownKeys, wparam);
            break;
        }
        case WM_CHAR:
            ::MultiByteToWideChar(CP_UTF8, MB_PRECOMPOSED, (char*)&wparam, 1, &wch, 1);
            if (KeyPressCallback) KeyPressCallback(wch);
        break;
        case WM_SIZE:
            windowWidth_  = LOWORD(lparam);
            windowHeight_ = HIWORD(lparam);
            UpdateRenderArea();
            if (WindowResizeCallback) WindowResizeCallback(windowWidth_, windowHeight_);
            break;
        case WM_MOVE:
            windowPosX_ = LOWORD(lparam);
            windowPosY_ = HIWORD(lparam);
            break;
        case WM_CLOSE:
        case WM_DESTROY:
            PostQuitMessage(0);
            break;
        default:
            result = DefWindowProcA(window, msg, wparam, lparam);
            break;
    }
    return result;
}

static HWND WindowCreate(HINSTANCE inst)
{
    WNDCLASSA window_class{};
    window_class.style         = CS_HREDRAW | CS_VREDRAW | CS_OWNDC;
    window_class.lpfnWndProc   = WindowCallback;
    window_class.hInstance     = inst;
    window_class.hCursor       = LoadCursor(0, IDC_ARROW);
    window_class.hbrBackground = 0;
    window_class.lpszClassName = "ASTLWindow";
    window_class.hIcon         = LoadIconA(inst, "duck_icon");
    
    if (!RegisterClassA(&window_class))
        FatalError("Failed to register window.");
    
    // Specify a desired width and height, then adjust the rect so the window's client area will be that size.
    RECT rect{};
    rect.right  = windowWidth_;
    rect.bottom = windowHeight_;
    const DWORD window_style = WS_OVERLAPPED     |
                               WS_CAPTION        |
                               WS_SYSMENU        |
                               WS_MINIMIZEBOX    |
                               WS_MAXIMIZEBOX;

    AdjustWindowRect(&rect, window_style, false);
    
    HWND window = CreateWindowExA(0,
                                  window_class.lpszClassName,
                                  WindowName, window_style,
                                  CW_USEDEFAULT, CW_USEDEFAULT,
                                  windowWidth_, windowHeight_,
                                  0, 0, inst, 0);
    
    if (!window) FatalError("Failed to create window.");
    return window;
}

float g_axDeltaTime = 0.001f;
float GetDeltaTime() { return g_axDeltaTime; }

int WINAPI WinMain(HINSTANCE inst, HINSTANCE prev, LPSTR cmd_line, int show)
{
    AXInit();
    
    hwnd        = WindowCreate(inst);
    HDC   dc    = GetDC(hwnd);
    HGLRC rc    = InitOpenGL(dc);
    
    gladLoaderLoadGL();
    
    ShowWindow(hwnd, show);
    UpdateWindow(hwnd);
    // first thing that we will see is going to be black color instead of white
    // if we clear before starting the engine
    glClearColor(0.3f, 0.3f, 0.3f, 1.0f); 
    SwapBuffers(dc);
    
    AXStart();

    LARGE_INTEGER frequency, prevTime, currentTime;
    QueryPerformanceFrequency(&frequency);
    QueryPerformanceCounter(&prevTime);

    bool running = true;
    while (running)
    {   
        MSG msg;
        while (PeekMessageA(&msg, 0, 0, 0, PM_REMOVE))
        {
            if (msg.message == WM_QUIT) 
            {
              running = false;
            }
            else 
            {
                TranslateMessage(&msg);
                DispatchMessageA(&msg);
            }
        }

        SetPressedAndReleasedKeys();
        
        QueryPerformanceCounter(&currentTime);
        g_axDeltaTime = (float)(currentTime.QuadPart - prevTime.QuadPart) / frequency.QuadPart;
        prevTime = currentTime;

        glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
        
        // Do OpenGL rendering here
        AXLoop();
        SwitchToThread();
        wglSwapIntervalEXT(1); // vsync on
        SwapBuffers(dc);

        RecordLastKeys();
    }

    AXExit();
    DestroyRenderer();
    wglMakeCurrent(dc, 0);
    ReleaseDC(hwnd, dc);
    wglDeleteContext(rc);
    DestroyWindow(hwnd);
    return 0;
}

#endif

void UpdateRenderArea()
{
#ifdef __ANDROID__
    EGLint width, height;
    eglQuerySurface(display_, surface_, EGL_WIDTH, &width);
    eglQuerySurface(display_, surface_, EGL_HEIGHT, &height);
    if (width != windowWidth_ || height != windowHeight_)
    {
      windowWidth_ = width;
      windowHeight_ = height;
    }
#else
    glViewport(0, 0, windowWidth_, windowHeight_);
#endif
}

void HandleInput()
{
#ifdef __ANDROID__
  // handle all queued inputs
  android_input_buffer* inputBuffer = android_app_swap_input_buffers(g_android_app);
  if (!inputBuffer) return; // no inputs yet.

  // handle motion events (motionEventsCounts can be 0).
  for (auto i = 0; i < inputBuffer->motionEventsCount; i++)
  {
    GameActivityMotionEvent& motionEvent = inputBuffer->motionEvents[i];
    int32_t action = motionEvent.action;

    // Find the pointer index, mask and bitshift to turn it into a readable value.
    int32_t pointerIndex = (action & AMOTION_EVENT_ACTION_POINTER_INDEX_MASK)
      >> AMOTION_EVENT_ACTION_POINTER_INDEX_SHIFT;
    AX_LOG("Pointer(s): ");

    // get the x and y position of this event if it is not ACTION_MOVE.
    GameActivityPointerAxes& pointer = motionEvent.pointers[pointerIndex];
    float x = GameActivityPointerAxes_getX(&pointer);
    float y = GameActivityPointerAxes_getY(&pointer);

    // determine the action type and process the event accordingly.
    switch (action & AMOTION_EVENT_ACTION_MASK)
    {
    case AMOTION_EVENT_ACTION_DOWN:
    case AMOTION_EVENT_ACTION_POINTER_DOWN:
      AX_LOG("( %i: %f, %f ) pointer down\n", pointer.id, x, y);
      break;

    case AMOTION_EVENT_ACTION_CANCEL:
      // treat the CANCEL as an UP event: doing nothing in the app, except
      // removing the pointer from the cache if pointers are locally saved.
      // code pass through on purpose.
    case AMOTION_EVENT_ACTION_UP:
    case AMOTION_EVENT_ACTION_POINTER_UP:
      AX_LOG("( %i: %f, %f ) pointer up\n", pointer.id, x, y);
      break;

    case AMOTION_EVENT_ACTION_MOVE:
      // There is no pointer index for ACTION_MOVE, only a snapshot of
      // all active pointers; app needs to cache previous active pointers
      // to figure out which ones are actually moved.
      for (auto index = 0; index < motionEvent.pointerCount; index++)
      {
        pointer = motionEvent.pointers[index];
        x = GameActivityPointerAxes_getX(&pointer);
        y = GameActivityPointerAxes_getY(&pointer);
        AX_LOG("( %i: %f, %f ) \n", pointer.id, x, y);

        if (index != (motionEvent.pointerCount - 1)) AX_LOG(",");
      }
      AX_LOG("Pointer Move");
      break;
    default:
      AX_LOG("Unknown MotionEvent Action: %i", action);
    }
    AX_LOG("Pointer Move\n");
  }
  // clear the motion input count in this buffer for main thread to re-use.
  android_app_clear_motion_events(inputBuffer);

  // handle input key events.
  for (auto i = 0; i < inputBuffer->keyEventsCount; i++)
  {
    auto& keyEvent = inputBuffer->keyEvents[i];
    AX_LOG("Key: %i ", keyEvent.keyCode);
    switch (keyEvent.action)
    {
    case AKEY_EVENT_ACTION_DOWN: AX_LOG("Key Down %i\n", keyEvent.action); break;
    case AKEY_EVENT_ACTION_UP:  AX_LOG("Key Up %i\n", keyEvent.action); break;
    case AKEY_EVENT_ACTION_MULTIPLE:
      // Deprecated since Android API level 29.
      AX_LOG("Multiple Key Actions %i\n", keyEvent.action);
      break;
    default:
      AX_LOG("Unknown KeyEvent Action: %i \n", keyEvent.action);
    }
  }
  // clear the key input count too.
  android_app_clear_key_events(inputBuffer);
#endif
}

void TerminateWindow()
{
#ifdef __ANDROID__
    eglMakeCurrent(display_, EGL_NO_SURFACE, EGL_NO_SURFACE, EGL_NO_CONTEXT);
    eglDestroyContext(display_, context_);
    eglDestroySurface(display_, surface_);
    eglTerminate(display_);
    display_ = EGL_NO_DISPLAY;
    surface_ = EGL_NO_SURFACE;
    context_ = EGL_NO_CONTEXT;
#endif
}