#pragma once
#include "Transform.hpp"
#include "../Window.hpp"

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

	Camera() { GetMonitorSize(&monitorSize.x, &monitorSize.y);  }

	Camera(Vector2i xviewPortSize) : viewportSize(xviewPortSize), 
	                                 position(MakeVec3(0.0f,4.0f, 15.0f)), 
	                                 targetPosition(MakeVec3(0.0f, 4.0f, 15.0f)), Front(MakeVec3(0.0f,0.0f,-1.0f))
	{
		GetMonitorSize(&monitorSize.x, &monitorSize.y);
		RecalculateProjection(xviewPortSize.x, xviewPortSize.y);
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

	void SetCursorPos(int x, int y) 
	{
		SetMousePos((float)x, (float)y); 
		mouseOld = MakeVec2((float)x, (float)y); 
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
		bool pressing = GetMouseDown(MouseButton_Right);
		if (!pressing) { wasPressing = false; return; }
	
		float dt = (float)GetDeltaTime() * 2.0f;
		float speed = dt * (1.0f + GetKeyDown(Key_SHIFT) * 2.0f) * 2.0f;
	
		Vector2f mousePos;
		GetMousePos(&mousePos.x, &mousePos.y);
		Vector2f diff = mousePos - mouseOld;
		
		if (wasPressing && diff.x + diff.y < 130.0f)
		{
			pitch -= diff.y * dt * senstivity;
			yaw   += diff.x * dt * senstivity;
			yaw   = FMod(yaw + 180.0f, 360.0f) - 180.0f;
			pitch = Clamp(pitch, -89.0f, 89.0f);
		}
		
		Front.x = Cos(yaw * DegToRad) * Cos(pitch * DegToRad);
		Front.y = Sin(pitch * DegToRad);
		Front.z = Sin(yaw * DegToRad) * Cos(pitch * DegToRad);
		Front.NormalizeSelf();
		// also re-calculate the Right and Up vector
		Right = Vector3f::Normalize(Vector3f::Cross(Front, Vector3f::Up()));  // normalize the vectors, because their length gets closer to 0 the more you look up or down which results in slower movement.
		Up = Vector3f::Normalize(Vector3f::Cross(Right, Front));
	
		if (GetKeyDown('D')) position += Right * speed;
		if (GetKeyDown('A')) position -= Right * speed;
		if (GetKeyDown('W')) position += Front * speed;
		if (GetKeyDown('S')) position -= Front * speed;
		if (GetKeyDown('Q')) position -= Up * speed;
		if (GetKeyDown('E')) position += Up * speed;
	
		mouseOld = mousePos;
		wasPressing = true;
	
		InfiniteMouse(mousePos);
		RecalculateView();
	}

	Ray ScreenPointToRay(Vector2f pos) const
	{
		Vector2f coord = MakeVec2(pos.x / (float)viewportSize.x, pos.y / (float)viewportSize.y);
		coord.y = 1.0f - coord.y;
		coord = coord * 2.0f - 1.0f;
		Vector4f target = Matrix4::Vector4Transform(MakeVec4(coord.x, coord.y, 1.0f, 1.0f), inverseProjection);
		target /= target.w;
		target = Matrix4::Vector4Transform(target, inverseView);
		Vector3f rayDir = Vector3f::Normalize(target.xyz());
		return MakeRay(position, rayDir);
	}

	RaySSE ScreenPointToRaySSE(Vector2f pos) const
	{
		Vector2f coord = MakeVec2(pos.x / (float)viewportSize.x, pos.y / (float)viewportSize.y);
		coord.y = 1.0f - coord.y;
		coord = coord * 2.0f - 1.0f;
		Vector4f target = Matrix4::Vector4Transform(MakeVec4(coord.x, coord.y, 1.0f, 1.0f), inverseProjection);
		target /= target.w;
		target = Matrix4::Vector4Transform(target, inverseView);
		Vector3f rayDir = Vector3f::Normalize(target.xyz());
		RaySSE ray;
		ray.origin = _mm_load_ps(&position.x);
		ray.direction = _mm_loadu_ps(&rayDir.x);
		ray.direction = _mm_insert_ps(ray.direction, _mm_setzero_ps(), 3);
		return ray;
	}
};
