
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
#include <stdlib.h> // exit failure

int windowWidth_   = 1240;
int windowHeight_  = 720;
int windowPosX_    = 500;
int windowPosY_    = 500;

static char WindowName[128]{ 'A', 'S', 'T', 'L' };
static void(*WindowMoveCallback)  (int, int) = nullptr;
static void(*WindowResizeCallback)(int, int) = nullptr;

extern void AXInit();
extern int  AXStart();
extern void AXLoop();
extern void AXExit();

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

void GetWindowSize(int* x, int* y) { }

void GetWindowPos(int* x, int* y) { }

void SetWindowName(const char* name) { }

#elif defined(AX_USE_WINDOW) && defined(_WIN32)

#define NOMINMAX
#define WIN32_LEAN_AND_MEAN 
#define VC_EXTRALEAN
#include <Windows.h>

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

void SetWindowResizeCallback(void(*callback)(int, int)) { WindowResizeCallback = callback; }
void SetWindowMoveCallback  (void(*callback)(int, int)) { WindowMoveCallback = callback; }
void GetWindowSize          (int* x, int* y)            { *x = windowWidth_; *y = windowHeight_; }
void GetWindowPos           (int* x, int* y)            { *x = windowPosX_; *y = windowPosY_; }

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

// See https://www.khronos.org/registry/OpenGL/extensions/ARB/WGL_ARB_create_context.txt for all values
// See https://www.khronos.org/registry/OpenGL/extensions/ARB/WGL_ARB_pixel_format.txt for all values
typedef HGLRC WINAPI wglCreateContextAttribsARB_type(HDC hdc, HGLRC hShareContext, const int* attribList);
wglCreateContextAttribsARB_type* wglCreateContextAttribsARB;

typedef BOOL WINAPI wglChoosePixelFormatARB_type(HDC hdc, const int* piAttribIList, const FLOAT* pfAttribFList, UINT nMaxFormats, int* piFormats, UINT* nNumFormats);
wglChoosePixelFormatARB_type* wglChoosePixelFormatARB;

static void fatal_error(char* msg)
{
  MessageBoxA(NULL, msg, "Error", MB_OK | MB_ICONEXCLAMATION);
  exit(EXIT_FAILURE);
}

static void init_opengl_extensions(void)
{
    // Before we can load extensions, we need a dummy OpenGL context, created using a dummy window.
    // We use a dummy window because you can only set the pixel format for a window once.
    WNDCLASSA window_class{};
    window_class.style         = CS_HREDRAW | CS_VREDRAW | CS_OWNDC;
    window_class.lpfnWndProc   = DefWindowProcA;
    window_class.hInstance     = GetModuleHandle(0);
    window_class.lpszClassName = "Dummy_WGL_djuasiodwa";
    
    if (!RegisterClass(&window_class)) 
        fatal_error("Failed to register dummy OpenGL window.");
    
    HWND dummy_window = CreateWindowExA(0, window_class.lpszClassName, "ASTL Window",
                                        0, CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT, 0,
                                        0, window_class.hInstance, 0);
    
    if (!dummy_window)
        fatal_error("Failed to create dummy OpenGL window.");
    
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
    if (!pixel_format) fatal_error("Failed to find a suitable pixel format.");
    
    if (!SetPixelFormat(dummy_dc, pixel_format, &pfd)) fatal_error("Failed to set the pixel format.");
    
    HGLRC dummy_context = wglCreateContext(dummy_dc);
    if (!dummy_context) fatal_error("Failed to create a dummy OpenGL rendering context.");
    
    if (!wglMakeCurrent(dummy_dc, dummy_context)) fatal_error("Failed to activate dummy OpenGL rendering context.");
    
    wglCreateContextAttribsARB = (wglCreateContextAttribsARB_type*)wglGetProcAddress("wglCreateContextAttribsARB");
    wglChoosePixelFormatARB    = (wglChoosePixelFormatARB_type*)wglGetProcAddress("wglChoosePixelFormatARB");
    
    wglMakeCurrent(dummy_dc, 0);
    wglDeleteContext(dummy_context);
    ReleaseDC(dummy_window, dummy_dc);
    DestroyWindow(dummy_window);
}

static HGLRC init_opengl(HDC real_dc)
{
    init_opengl_extensions();
    // Now we can choose a pixel format the modern way, using wglChoosePixelFormatARB.
    int pixel_format_attribs[] = {
        0x2001,     GL_TRUE, // WGL_DRAW_TO_WINDOW_ARB
        0x2010,     GL_TRUE, // WGL_SUPPORT_OPENGL_ARB
        0x2011,     GL_TRUE, // WGL_DOUBLE_BUFFER_ARB
        0x2003,     0x2027,  // WGL_ACCELERATION_ARB, WGL_FULL_ACCELERATION_ARB
        0x2013,     0x202B,  // WGL_PIXEL_TYPE_ARB, WGL_TYPE_RGBA_ARB
        0x2014,         32,  // WGL_COLOR_BITS_ARB
        0x2022,         24,  // WGL_DEPTH_BITS_ARB
        0x2023,          8,  // WGL_STENCIL_BITS_ARB
        0
    };
    
    int pixel_format;
    UINT num_formats;
    wglChoosePixelFormatARB(real_dc, pixel_format_attribs, 0, 1, &pixel_format, &num_formats);
    if (!num_formats) 
      fatal_error("Failed to set the OpenGL 3.3 pixel format.");
    
    PIXELFORMATDESCRIPTOR pfd;
    DescribePixelFormat(real_dc, pixel_format, sizeof(pfd), &pfd);
    if (!SetPixelFormat(real_dc, pixel_format, &pfd)) 
      fatal_error("Failed to set the OpenGL 3.3 pixel format.");
    
    // Specify that we want to create an OpenGL 3.2 core profile context
    int gl32_attribs[] = {
        0x2091, 3, // WGL_CONTEXT_MAJOR_VERSION_ARB
        0x2092, 2, // WGL_CONTEXT_MINOR_VERSION_ARB
        0x9126,  0x00000001, // WGL_CONTEXT_PROFILE_MASK_ARB, WGL_CONTEXT_CORE_PROFILE_BIT_ARB
        0,
    };
    
    HGLRC gl32_context = wglCreateContextAttribsARB(real_dc, 0, gl32_attribs);
    if (!gl32_context) 
      fatal_error("Failed to create OpenGL 3.2 context.");
    
    if (!wglMakeCurrent(real_dc, gl32_context)) 
      fatal_error("Failed to activate OpenGL 3.2 rendering context.");
    
    return gl32_context;
}

static LRESULT CALLBACK window_callback(HWND window, UINT msg, WPARAM wparam, LPARAM lparam)
{
    LRESULT result = 0;
    
    switch (msg)
    {
    case WM_SIZE:
        windowWidth_  = LOWORD(lparam);
        windowHeight_ = HIWORD(lparam);
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

static HWND create_window(HINSTANCE inst)
{
    WNDCLASSA window_class{};
    window_class.style         = CS_HREDRAW | CS_VREDRAW | CS_OWNDC;
    window_class.lpfnWndProc   = window_callback;
    window_class.hInstance     = inst;
    window_class.hCursor       = LoadCursor(0, IDC_ARROW);
    window_class.hbrBackground = 0;
    window_class.lpszClassName = "ASTLWindow";
    
    if (!RegisterClassA(&window_class)) fatal_error("Failed to register window.");
    
    // Specify a desired width and height, then adjust the rect so the window's client area will be
    // that size.
    RECT rect{};
    rect.right  = windowWidth_;
    rect.bottom = windowHeight_;
    DWORD window_style = WS_OVERLAPPEDWINDOW;
    AdjustWindowRect(&rect, window_style, false);
    
    HWND window = CreateWindowExA(
      0,
      window_class.lpszClassName,
      WindowName,
      window_style,
      CW_USEDEFAULT,
      CW_USEDEFAULT,
      windowWidth_,
      windowHeight_,
      0,
      0,
      inst,
      0);
    
    if (!window) fatal_error("Failed to create window.");
    return window;
}

int WINAPI WinMain(HINSTANCE inst, HINSTANCE prev, LPSTR cmd_line, int show)
{
    AXInit();

    hwnd        = create_window(inst);
    HDC   dc    = GetDC(hwnd);
    HGLRC rc    = init_opengl(dc);
    
    gladLoaderLoadGL();
    ShowWindow(hwnd, show);
    UpdateWindow(hwnd);
    AXStart();
    
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
        
        glClearColor(1.0f, 0.5f, 0.5f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        // Do OpenGL rendering here
        AXLoop();
        SwapBuffers(dc);
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