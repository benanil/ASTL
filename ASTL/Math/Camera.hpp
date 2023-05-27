#pragma once
#include "Transform.hpp"
#include "Vector4.hpp"
#include "../Window.hpp"
#include "../Logger.hpp"

struct Camera
{
	Matrix4 projection;
	Matrix4 view;
	
	Matrix4 inverseProjection;
	Matrix4 inverseView;

	float verticalFOV = 65.0f;
	float nearClip = 0.01f;
	float farClip = 500.0f;
	
	Vector2i viewportSize, monitorSize;
	
	Vector3f position, targetPosition;
	Vector2f mouseOld;
	
	Vector3f Front, Right, Up;
	
	float pitch = 0.0f, yaw = -90.0f , senstivity = 20.0f;
	int projWidth, projHeight;

	bool wasPressing = false;

	Camera() {}

	Camera(Vector2i xviewPortSize)
	: viewportSize(xviewPortSize), position(0.0f,4.0f, 15.0f), targetPosition(0.0f, 4.0f, 15.0f), Front(0.0f,0.0f,-1.0f)
	{
		monitorSize = Window::GetMonitorScale();
		RecalculateProjection(xviewPortSize.x, xviewPortSize.y);
		RecalculateView();
	}
	
	void SetCursorPos(int x, int y) 
	{
		Window::SetMouseScreenPos(Vector2i(x, y)); 
		mouseOld = Vector2f((float)x, (float)y); 
	}

	void InfiniteMouse(const Vector2f& point)
	{
		if (point.x > monitorSize.x - 2) SetCursorPos(3, (int)point.y);
		if (point.y > monitorSize.y - 2) SetCursorPos((int)point.x, 3);
		
		if (point.x < 2) SetCursorPos(monitorSize.x - 3, (int)point.y);
		if (point.y < 2) SetCursorPos((int)point.x, monitorSize.y - 3);
	}

	void Update()
	{
		bool pressing = Window::GetMouseButton(MouseButton_Right);
		if (!pressing) { wasPressing = false; return; }

		float dt = (float)Window::DeltaTime();
		float speed = dt * (1.0f + Window::GetKey(KeyCode_LEFT_SHIFT) * 2.0f) * 2.0f;

		const Vector2f mousePos = ToVector2f(Window::GetMouseScreenPos());
		Vector2f diff = mousePos - mouseOld;
		
		if (wasPressing && diff.x + diff.y < 130.0f)
		{
			pitch -= diff.y * dt * senstivity;
			yaw   += diff.x * dt * senstivity;
			pitch = Clamp(pitch, -89.0f, 89.0f);
		}
		
		Front.x = Cos(yaw * DegToRad) * Cos(pitch * DegToRad);
		Front.y = Sin(pitch * DegToRad);
		Front.z = Sin(yaw * DegToRad) * Cos(pitch * DegToRad);
		Front.NormalizeSelf();
		// also re-calculate the Right and Up vector
		Right = Vector3f::Normalize(Vector3f::Cross(Front, Vector3f::Up()));  // normalize the vectors, because their length gets closer to 0 the more you look up or down which results in slower movement.
		Up = Vector3f::Normalize(Vector3f::Cross(Right, Front));

		if (Window::GetKey(KeyCode_D)) position += Right * speed;
		if (Window::GetKey(KeyCode_A)) position -= Right * speed;
		if (Window::GetKey(KeyCode_W)) position += Front * speed;
		if (Window::GetKey(KeyCode_S)) position -= Front * speed;
		if (Window::GetKey(KeyCode_Q)) position -= Up * speed;
		if (Window::GetKey(KeyCode_E)) position += Up * speed;

		mouseOld = mousePos;
		wasPressing = true;

		InfiniteMouse(mousePos);
		RecalculateView();
	}

	void RecalculateProjection(int width, int height)
	{
		projWidth = width; projHeight = height;
		projection = Matrix4::PerspectiveFovRH(verticalFOV * DegToRad, (float)width, (float)height, nearClip, farClip);
		inverseProjection = Matrix4::Inverse(projection);	
	}

	void RecalculateView()
	{
		view = Matrix4::LookAtRH(position, Front, Vector3f::Up());
		inverseView = Matrix4::Inverse(view);
	}

	Ray ScreenPointToRay(Vector2f pos) const
	{
		Vector2f coord(pos.x / (float)viewportSize.x, pos.y / (float)viewportSize.y);
		coord.y = 1.0f - coord.y;
		coord = coord * 2.0f - 1.0f;
		Vector4 target = Matrix4::Vector4Transform(Vector4(coord.x, coord.y, 1.0f, 1.0f), inverseProjection);
		target /= target.w;
		target = Matrix4::Vector4Transform(target, inverseView);
		Vector3f rayDir = Vector3f::Normalize(target.xyz());
		return Ray(position, rayDir);
	}

	RaySSE ScreenPointToRaySSE(Vector2f pos) const
	{
		Vector2f coord(pos.x / (float)viewportSize.x, pos.y / (float)viewportSize.y);
		coord.y = 1.0f - coord.y;
		coord = coord * 2.0f - 1.0f;
		Vector4 target = Matrix4::Vector4Transform(Vector4(coord.x, coord.y, 1.0f, 1.0f), inverseProjection);
		target /= target.w;
		target = Matrix4::Vector4Transform(target, inverseView);
		Vector3f rayDir = Vector3f::Normalize(target.xyz());
		RaySSE ray;
		ray.origin = _mm_load_ps(&position.x);
		ray.direction = _mm_loadu_ps(&rayDir.x);
		ray.direction.m128_f32[3] = 0;
		return ray;
	}
};
