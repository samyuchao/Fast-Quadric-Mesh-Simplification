//------------------------------------------------------------------------------
// <auto-generated />
//
// This file was automatically generated by SWIG (http://www.swig.org).
// Version 3.0.12
//
// Do not make changes to this file unless you know what you are doing--modify
// the SWIG interface file instead.
//------------------------------------------------------------------------------


public class vec3f : global::System.IDisposable {
  private global::System.Runtime.InteropServices.HandleRef swigCPtr;
  protected bool swigCMemOwn;

  internal vec3f(global::System.IntPtr cPtr, bool cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = new global::System.Runtime.InteropServices.HandleRef(this, cPtr);
  }

  internal static global::System.Runtime.InteropServices.HandleRef getCPtr(vec3f obj) {
    return (obj == null) ? new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero) : obj.swigCPtr;
  }

  ~vec3f() {
    Dispose();
  }

  public virtual void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != global::System.IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          cppdllPINVOKE.delete_vec3f(swigCPtr);
        }
        swigCPtr = new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero);
      }
      global::System.GC.SuppressFinalize(this);
    }
  }

  public double x {
    set {
      cppdllPINVOKE.vec3f_x_set(swigCPtr, value);
    } 
    get {
      double ret = cppdllPINVOKE.vec3f_x_get(swigCPtr);
      return ret;
    } 
  }

  public double y {
    set {
      cppdllPINVOKE.vec3f_y_set(swigCPtr, value);
    } 
    get {
      double ret = cppdllPINVOKE.vec3f_y_get(swigCPtr);
      return ret;
    } 
  }

  public double z {
    set {
      cppdllPINVOKE.vec3f_z_set(swigCPtr, value);
    } 
    get {
      double ret = cppdllPINVOKE.vec3f_z_get(swigCPtr);
      return ret;
    } 
  }

  public vec3f() : this(cppdllPINVOKE.new_vec3f__SWIG_0(), true) {
  }

  public vec3f(vector3 a) : this(cppdllPINVOKE.new_vec3f__SWIG_1(vector3.getCPtr(a)), true) {
    if (cppdllPINVOKE.SWIGPendingException.Pending) throw cppdllPINVOKE.SWIGPendingException.Retrieve();
  }

  public vec3f(double X, double Y, double Z) : this(cppdllPINVOKE.new_vec3f__SWIG_2(X, Y, Z), true) {
  }

  public vec3f v3() {
    vec3f ret = new vec3f(cppdllPINVOKE.vec3f_v3(swigCPtr), true);
    return ret;
  }

  public double dot(vec3f a) {
    double ret = cppdllPINVOKE.vec3f_dot(swigCPtr, vec3f.getCPtr(a));
    if (cppdllPINVOKE.SWIGPendingException.Pending) throw cppdllPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public vec3f cross(vec3f a, vec3f b) {
    vec3f ret = new vec3f(cppdllPINVOKE.vec3f_cross(swigCPtr, vec3f.getCPtr(a), vec3f.getCPtr(b)), true);
    if (cppdllPINVOKE.SWIGPendingException.Pending) throw cppdllPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public double angle(vec3f v) {
    double ret = cppdllPINVOKE.vec3f_angle(swigCPtr, vec3f.getCPtr(v));
    if (cppdllPINVOKE.SWIGPendingException.Pending) throw cppdllPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public double angle2(vec3f v, vec3f w) {
    double ret = cppdllPINVOKE.vec3f_angle2(swigCPtr, vec3f.getCPtr(v), vec3f.getCPtr(w));
    if (cppdllPINVOKE.SWIGPendingException.Pending) throw cppdllPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public vec3f rot_x(double a) {
    vec3f ret = new vec3f(cppdllPINVOKE.vec3f_rot_x(swigCPtr, a), true);
    return ret;
  }

  public vec3f rot_y(double a) {
    vec3f ret = new vec3f(cppdllPINVOKE.vec3f_rot_y(swigCPtr, a), true);
    return ret;
  }

  public void clamp(double min, double max) {
    cppdllPINVOKE.vec3f_clamp(swigCPtr, min, max);
  }

  public vec3f rot_z(double a) {
    vec3f ret = new vec3f(cppdllPINVOKE.vec3f_rot_z(swigCPtr, a), true);
    return ret;
  }

  public vec3f invert() {
    vec3f ret = new vec3f(cppdllPINVOKE.vec3f_invert(swigCPtr), true);
    return ret;
  }

  public vec3f frac() {
    vec3f ret = new vec3f(cppdllPINVOKE.vec3f_frac(swigCPtr), true);
    return ret;
  }

  public vec3f integer() {
    vec3f ret = new vec3f(cppdllPINVOKE.vec3f_integer(swigCPtr), true);
    return ret;
  }

  public double length() {
    double ret = cppdllPINVOKE.vec3f_length(swigCPtr);
    return ret;
  }

  public vec3f normalize(double desired_length) {
    vec3f ret = new vec3f(cppdllPINVOKE.vec3f_normalize__SWIG_0(swigCPtr, desired_length), true);
    return ret;
  }

  public vec3f normalize() {
    vec3f ret = new vec3f(cppdllPINVOKE.vec3f_normalize__SWIG_1(swigCPtr), true);
    return ret;
  }

}
