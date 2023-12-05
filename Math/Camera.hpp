#pragma once
#include "Transform.hpp"

struct Camera
{
	Matrix4 projection;
	Matrix4 view;
	
	Matrix4 inverseProjection;
	Matrix4 inverseView;

	float verticalFOV = 65.0f;
	float nearClip = 0.01f;
	float farClip = 500.0f;
	
	Vector2i viewportSize;
	
	Vector3f position{0.0f,4.0f, 15.0f};
	Vector3f Front{0.0f,0.0f,-1.0f};

	Camera() {}

	Camera(Vector2i xviewPortSize) : viewportSize(xviewPortSize)
	{
		RecalculateProjection(xviewPortSize.x, xviewPortSize.y);
		RecalculateView();
	}
	
	void RecalculateProjection(int width, int height)
	{
		viewportSize.x = width; viewportSize.y = height;
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
		Vector2f coord{pos.x / (float)viewportSize.x, pos.y / (float)viewportSize.y};
		coord.y = 1.0f - coord.y;
		coord = coord * 2.0f - 1.0f;
		Vector4f target = Matrix4::Vector4Transform(MakeVec4(coord.x, coord.y, 1.0f, 1.0f), inverseProjection);
		target /= target.w;
		target = Matrix4::Vector4Transform(target, inverseView);
		Vector3f rayDir = Vector3f::Normalize(target.xyz());
		return MakeRay(position, rayDir);
	}

#if defined(AX_SUPPORT_SSE)
	RaySSE ScreenPointToRaySSE(Vector2f pos) const
	{
		Vector2f coord(pos.x / (float)viewportSize.x, pos.y / (float)viewportSize.y);
		coord.y = 1.0f - coord.y;
		coord = coord * 2.0f - 1.0f;
		Vector4f target = Matrix4::Vector4Transform(Vector4f(coord.x, coord.y, 1.0f, 1.0f), inverseProjection);
		target /= target.w;
		target = Matrix4::Vector4Transform(target, inverseView);
		Vector3f rayDir = Vector3f::Normalize(target.xyz());
		RaySSE ray;
		ray.origin = _mm_load_ps(&position.x);
		ray.direction = _mm_loadu_ps(&rayDir.x);
		ray.direction.m128_f32[3] = 0;
		return ray;
	}
#endif
};
