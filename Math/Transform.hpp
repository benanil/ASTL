
#pragma once

#include "Matrix.hpp"

struct Transform
{
public:
	Matrix4 transform = Matrix4::Identity();
	Quaternion rotation{};
	Vector3f scale {1, 1, 1};
	bool needsUpdate = false;

	Transform() {}

	void SetScale(Vector3f scale) {
		this->scale = scale; needsUpdate = true;
	}
	
	void SetPosition(float x, float y, float z)
	{
		transform.m[3][0] = x;
		transform.m[3][1] = y;
		transform.m[3][2] = z;
	}

	void SetPosition(Vector3f position) 
	{
		SetPosition(position.x, position.y, position.z);
	}
	
	Vector3f GetPosition() const
	{
		return { transform.m[3][0], transform.m[3][1], transform.m[3][2] };
	}

	void SetRotationEuler(Vector3f euler) 
	{
		this->rotation = Quaternion::FromEuler(euler); needsUpdate = true;
	}
	
	void SetRotationEulerDegree(Vector3f euler) 
	{
		this->rotation = Quaternion::FromEuler(euler * DegToRad); needsUpdate = true;
	}

	void SetRotationQuaternion(const Quaternion& rotation) 
	{
		this->rotation = rotation; needsUpdate = true;
	}

	void SetMatrix(const Matrix4& matrix)
	{
		this->rotation = Matrix4::ExtractRotation(matrix);
		this->scale    = Matrix4::ExtractScale(matrix);
		this->transform = matrix;
	}

	void CalculateMatrix()
	{
		transform = Matrix4::FromPosition(GetPosition()) * Matrix4::FromQuaternion(rotation) * Matrix4::CreateScale(scale);
	}

	Matrix4& GetMatrix()
	{
		if (needsUpdate) {
			transform = Matrix4::FromPosition(GetPosition()) * Matrix4::FromQuaternion(rotation) * Matrix4::CreateScale(scale);
		}
		return transform;
	}

	Quaternion GetRotation() const { return rotation; }
	Vector3f GetEulerDegree() const { return Quaternion::ToEulerAngles(rotation) * RadToDeg; }
	Vector3f GetEuler() const { return Quaternion::ToEulerAngles(rotation); }
	
	Vector3f GetForward() const { return rotation.GetForward(); }
	Vector3f GetUp()      const { return rotation.GetUp(); }
	Vector3f GetRight()   const { return rotation.GetRight(); }
};