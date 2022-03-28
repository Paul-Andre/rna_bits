struct V3 {
  float x,y,z;
};

V3 operator+(const V3 &a, const V3 &b) {
  V3 ret;
  ret.x = a.x + b.x;
  ret.y = a.y + b.y;
  ret.z = a.z + b.z;
  return ret;
}

V3 operator*(float a, const V3 &b) {
  V3 ret;
  ret.x = a * b.x;
  ret.y = a * b.y;
  ret.z = a * b.z;
  return ret;
}

V3 operator-(const V3 &a, const V3 &b) {
  V3 ret;
  ret.x = a.x - b.x;
  ret.y = a.y - b.y;
  ret.z = a.z - b.z;
  return ret;
}

V3 operator(const V3 &a, const V3 &b) {
  V3 ret;
  ret.x = a.x - b.x;
  ret.y = a.y - b.y;
  ret.z = a.z - b.z;
  return ret;
}







