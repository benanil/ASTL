
#pragma once
#include "Matrix.hpp"

class Transform
{
private:
	Matrix4 transform = Matrix4::Identity();
public:
	Vector3f position{};
	bool needsUpdate = false;
	Quaternion rotation{};
	Vector3f scale {1, 1, 1};

	Transform() {}

	void SetScale(Vector3f scale) {
		this->scale = scale; needsUpdate = true;
	}
	
	void SetPosition(Vector3f position) { this->position = position; needsUpdate = true; }
	
	void SetPosition(float x, float y, float z) { this->position.x = x; this->position.y = y; this->position.z = z; needsUpdate = true; }

	void SetRotationEuler(Vector3f euler) {
		this->rotation = Quaternion::FromEuler(euler); needsUpdate = true;
	}
	
	void SetRotationEulerDegree(Vector3f euler) {
		this->rotation = Quaternion::FromEuler(euler * DegToRad); needsUpdate = true;
	}

	void SetRotationQuaternion(const Quaternion& rotation) {
		this->rotation = rotation; needsUpdate = true;
	}

	void SetMatrix(const Matrix4& matrix)
	{
		this->rotation = Matrix4::ExtractRotation(matrix);
		this->position = Matrix4::ExtractPosition(matrix);
		this->scale    = Matrix4::ExtractScale(matrix);
		this->transform = matrix;
	}

	void UpdateMatrix()
	{
		transform = Matrix4::Identity() * Matrix4::FromPosition(position) * Matrix4::FromQuaternion(rotation) * Matrix4::CreateScale(scale);
	}

	void UpdatePosition()
	{
		transform.m[3][0] = position.x;
		transform.m[3][1] = position.y;
		transform.m[3][2] = position.z;
	}

	Matrix4& GetMatrix()
	{
		if (needsUpdate) {
			transform = Matrix4::Identity() * Matrix4::FromPosition(position) * Matrix4::FromQuaternion(rotation) * Matrix4::CreateScale(scale);
		}
		return transform;
	}

	Quaternion GetRotation() const { return rotation; }
	Vector3f GetEulerDegree() const { return Quaternion::ToEulerAngles(rotation) * RadToDeg; }
	Vector3f GetForward() const { return Vector3f(transform.m[2][0], transform.m[2][1], transform.m[2][2]); }
	Vector3f GetUp()      const { return Vector3f(transform.m[1][0], transform.m[1][1], transform.m[1][2]); }
	Vector3f GetRight()   const { return Vector3f(transform.m[0][0], transform.m[0][1], transform.m[0][2]); }
};