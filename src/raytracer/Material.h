#pragma once
//Contains the PBRmaterial class that holds information of the owning primitive 
//required for the sampling and attenuation of the hitting ray (e.g. albedo,
//roughness, emissiveness etc. either as scalars or maps if available), as well
//as the texture class
//PBR = physically based rendering

#include <memory>
#include <optional>
#include <vector>

#include "UtilConcepts.h"
#include "Vector.h"
#include "Image.h"
#include "Primitives.h"

template<typename T, Float_t space_type, Float_t color_type>
class Texture {
public:
  std::shared_ptr<Image<T>> im;

  Texture() = delete;
  Texture(std::shared_ptr<Image<T>> im) : im(im) {};
};

template<Float_t space_type>
struct UvMapTri {
  Vec2<space_type> uv_verts[3];
};
template<Float_t space_type>
struct UvMapSphere {
  Vec3<space_type> e_y;
};

enum Emissiveness_type {COSINE, COSINE_SQUARE, DISCRETE_BEAM, UNIFORM};

//Contains information to determine PBR properties of ascociated primitive
template<Float_t space_type, Float_t color_type>
class PBRMaterial {
public:
  //General pattern is for each propertie there is a constant value with a
  //default and a shrd ptr to a texture with an the ascociated UV map for
  //In case the shrd ptr of the specific texture is == nullptr the constant
  //value is used otherwise the texture will be used alongisde either of the
  //required Uv maps (which depends on the used  priomitive type e.g. tri or
  //sphere)
  //The Uv map is assumed to be set if it is required for Uv mapping. The 
  //program will crash if it isnt
  
  Pixel<color_type> albedo_const = Pixel<color_type>(0.5,0.5,0.5);
  std::shared_ptr<Texture<Pixel<color_type>,space_type,color_type>>
      albedo_map = nullptr;
  std::optional<UvMapTri<space_type>> albedo_uv_map_tri;
  std::optional<UvMapSphere<space_type>> albedo_uv_map_sphere;

  space_type roughness_const = 0.2;
  std::shared_ptr<Texture<space_type,space_type,color_type>>
      roughness_map = nullptr;
  std::optional<UvMapTri<space_type>> roughness_uv_map_tri;
  std::optional<UvMapSphere<space_type>> roughness_uv_map_sphere;

  space_type metalness_const = 0.0;
  std::shared_ptr<Texture<space_type,space_type,color_type>>
      metalness_map = nullptr;
  std::optional<UvMapTri<space_type>> metalness_uv_map_tri;
  std::optional<UvMapSphere<space_type>> metalness_uv_map_sphere;

  //The default normal is the geometry normal
  std::shared_ptr<Texture<Vec3<space_type>,space_type,color_type>>
      normal_map = nullptr;
  std::optional<UvMapTri<space_type>> normal_uv_map_tri;
  std::optional<UvMapSphere<space_type>> normal_uv_map_sphere;
  
  //For type of ditribution of emitted light
  Emissiveness_type emit_type = COSINE;
  std::optional<space_type> beam_width;
  //For color and base intensity
  Pixel<color_type> emitance_const = Pixel<color_type>(0.0,0.0,0.0);
  std::shared_ptr<Texture<Pixel<color_type>,space_type,color_type>>
      emitance_map = nullptr;
  std::optional<UvMapTri<space_type>> emitance_uv_map_tri;
  std::optional<UvMapSphere<space_type>> emitance_uv_map_sphere;

  //Translucency is currently handled by getting the transmissty value of the 
  //ray intersection point and then with a likely hood = to that refract the ray
  //based on snells law and otherwise generate the ray based on the normal alg
  space_type refractive_index_const = 1.0;
  std::shared_ptr<Texture<space_type,space_type,color_type>>
      refractive_index_map = nullptr;
  std::optional<UvMapTri<space_type>> refractive_index_uv_map_tri;
  std::optional<UvMapSphere<space_type>> refractive_index_uv_map_sphere;

  space_type transmissity_const = 0.0;
  std::shared_ptr<Texture<space_type,space_type,color_type>>
      transmissity_map = nullptr;
  std::optional<UvMapTri<space_type>> transmissity_uv_map_tri;
  std::optional<UvMapSphere<space_type>> transmissity_uv_map_sphere;
  
  PBRMaterial() {};
  PBRMaterial(
      const Pixel<color_type>& albedo,
      space_type roughness,
      space_type metalness,
      Emissiveness_type emit_type,
      std::optional<space_type> beam_width,
      const Pixel<color_type>& emitance,
      space_type refractive_index,
      space_type transmissity) : albedo_const(albedo),
                                 roughness_const(roughness),
                                 metalness_const(metalness),
                                 emit_type(emit_type),
                                 beam_width(beam_width),
                                 emitance_const(emitance),
                                 refractive_index_const(refractive_index),
                                 transmissity_const(transmissity) {};
  PBRMaterial(
      const PBRMaterial<space_type,color_type>& other) :
                                 albedo_const(other.albedo_const),
                                 roughness_const(other.roughness_const),
                                 metalness_const(other.metalness_const),
                                 emit_type(other.emit_type),
                                 beam_width(beam_width),
                                 emitance_const(other.emitance_const),
                                 refractive_index_const(other.refractive_index_const),
                                 transmissity_const(other.transmissity_const) {
    albedo_map = other.albedo_map;
    if (other.albedo_uv_map_tri) {
      albedo_uv_map_tri = other.albedo_uv_map_tri;
    }
    if (other.albedo_uv_map_sphere) {
      albedo_uv_map_sphere = other.albedo_uv_map_sphere;
    }

    roughness_map = other.roughness_map;
    if (other.roughness_uv_map_tri) {
      roughness_uv_map_tri = other.roughness_uv_map_tri;
    }
    if (other.roughness_uv_map_sphere) {
      roughness_uv_map_sphere = other.roughness_uv_map_sphere;
    }

    metalness_map = other.metalness_map;
    if (other.metalness_uv_map_tri) {
      metalness_uv_map_tri = other.metalness_uv_map_tri;
    }
    if (other.metalness_uv_map_sphere) {
      metalness_uv_map_sphere = other.metalness_uv_map_sphere;
    }

    normal_map = other.normal_map;
    if (other.normal_uv_map_tri) {
      normal_uv_map_tri = other.normal_uv_map_tri;
    }
    if (other.normal_uv_map_sphere) {
      normal_uv_map_sphere = other.normal_uv_map_sphere;
    }

    emitance_map = other.emitance_map;
    if (other.emitance_uv_map_tri) {
      emitance_uv_map_tri = other.emitance_uv_map_tri;
    }
    if (other.emitance_uv_map_sphere) {
      emitance_uv_map_sphere = other.emitance_uv_map_sphere;
    }

    refractive_index_map = other.refractive_index_map;
    if (other.refractive_index_uv_map_tri) {
      refractive_index_uv_map_tri = other.refractive_index_uv_map_tri;
    }
    if (other.refractive_index_uv_map_sphere) {
      refractive_index_uv_map_sphere = other.refractive_index_uv_map_sphere;
    }

    transmissity_map = other.transmissity_map;
    if (other.transmissity_uv_map_tri) {
      transmissity_uv_map_tri = other.transmissity_uv_map_tri;
    }
    if (other.transmissity_uv_map_sphere) {
      transmissity_uv_map_sphere = other.transmissity_uv_map_sphere;
    }
  }

  Pixel<color_type> GetAlbedoTri(
      const Vec3<space_type>& point_on_tri,
      const Tri<space_type,color_type>& tri) {
    if (!albedo_uv_map_tri) {
      return albedo_const;
    }

    return ExtractPixelFromTextureMap(
        albedo_uv_map_tri.value(),
        tri,
        point_on_tri,
        *(albedo_map.get()));
  }
  Pixel<color_type> GetAlbedoSphere(
      const Vec3<space_type>& point_on_surface,
      const Vec3<space_type>& sphere_origin) {
    if (!albedo_uv_map_sphere) {
      return albedo_const;
    }

    return ExtractPixelFromTextureMap(
        point_on_surface,
        sphere_origin,
        albedo_uv_map_sphere.value().e_y,
        *(albedo_map.get()));
  }

  space_type GetRoughnessTri(
      const Vec3<space_type>& point_on_tri,
      const Tri<space_type,color_type>& tri) {
    if (!roughness_uv_map_tri) {
      return roughness_const;
    }

    return ExtractScalarFromTextureMap(
        roughness_uv_map_tri.value(),
        tri,
        point_on_tri,
        *(roughness_map.get()));
  }
  space_type GetRoughnessSphere(
      const Vec3<space_type>& point_on_surface,
      const Vec3<space_type>& sphere_origin) {
    if (!roughness_uv_map_sphere) {
      return roughness_const;
    }

    return ExtractScalarFromTextureMap(
        point_on_surface,
        sphere_origin,
        roughness_uv_map_sphere.value().e_y,
        *(roughness_map.get()));
  }

  space_type GetMetalnessTri(
      const Vec3<space_type>& point_on_tri,
      const Tri<space_type,color_type>& tri) {
    if (!metalness_uv_map_tri) {
      return metalness_const;
    }

    return ExtractScalarFromTextureMap(
        metalness_uv_map_tri.value(),
        tri,
        point_on_tri,
        *(metalness_map.get()));
  }
  space_type GetMetalnessSphere(
      const Vec3<space_type>& point_on_surface,
      const Vec3<space_type>& sphere_origin) {
    if (!metalness_uv_map_sphere) {
      return metalness_const;
    }

    return ExtractScalarFromTextureMap(
        point_on_surface,
        sphere_origin,
        metalness_uv_map_sphere.value().e_y,
        *(metalness_map.get()));
  }

  Vec3<space_type> GetNormalTri(
      const Vec3<space_type>& point_on_tri,
      const Tri<space_type,color_type>& tri) {
    if (!normal_uv_map_tri) {
      return tri.geometry_normal;
    }

    return ExtractVec3FromTextureMap(
        normal_uv_map_tri.value(),
        tri,
        point_on_tri,
        *(normal_map.get()));
  }
  Vec3<space_type> GetNormalSphere(
      const Vec3<space_type>& point_on_surface,
      const Vec3<space_type>& sphere_origin) {
    if (!normal_uv_map_sphere) {
      Vec3<space_type> result{(point_on_surface - sphere_origin)};
      result.Normalize();
      return result;
    }

    return ExtractVec3FromTextureMap(
        point_on_surface,
        sphere_origin,
        normal_uv_map_sphere.value().e_y,
        *(normal_map.get()));
  }

  space_type GetEmmissiveFactor(space_type angle) {
    switch (emit_type) {
      case COSINE:
        {
        return std::cos(angle);
        }
      case COSINE_SQUARE:
        {
        space_type tmp = std::cos(angle);
        return tmp*tmp;
        }
      case DISCRETE_BEAM:
        {
        return angle <= beam_width ? 1.0 : 0;
        }
      case UNIFORM:
        {
        return 1.0;
        }
      default:
        {
        std::cout << "Unsupported emissive type!!!!" << std::endl;
        std::terminate();
        }
    }
    return 0.0; //unreachable
  }

  Pixel<color_type> GetEmitanceTri(
      const Vec3<space_type>& point_on_tri,
      const Tri<space_type,color_type>& tri,
      space_type ray_angle) {
    Pixel<color_type> result;
    if (!emitance_uv_map_tri) {
      result = emitance_const;
    } else {
      result = ExtractPixelFromTextureMap(
          emitance_uv_map_tri.value(),
          tri,
          point_on_tri,
          *(emitance_map.get()));
    }

    return result * GetEmmissiveFactor(ray_angle);
  }

  Pixel<color_type> GetEmitanceSphere(
      const Vec3<space_type>& point_on_surface,
      const Vec3<space_type>& sphere_origin,
      space_type ray_angle) {
    Pixel<color_type> result;
    if (!emitance_uv_map_sphere) {
      result = emitance_const;
    } else {
      result = ExtractPixelFromTextureMap(
        point_on_surface,
        sphere_origin,
        emitance_uv_map_sphere.value().e_y,
        *(emitance_map.get()));
    }

    return result * GetEmmissiveFactor(ray_angle);
  }

  space_type GetRefractiveIndexTri(
      const Vec3<space_type>& point_on_tri,
      const Tri<space_type,color_type>& tri) {
    if (!refractive_index_uv_map_tri) {
      return refractive_index_const;
    }

    return ExtractScalarFromTextureMap(
        refractive_index_uv_map_tri.value(),
        tri,
        point_on_tri,
        *(refractive_index_map.get()));
  }
  space_type GetRefractiveIndexSphere(
      const Vec3<space_type>& point_on_surface,
      const Vec3<space_type>& sphere_origin) {
    if (!refractive_index_uv_map_sphere) {
      return refractive_index_const;
    }

    return ExtractScalarFromTextureMap(
        point_on_surface,
        sphere_origin,
        refractive_index_uv_map_sphere.value().e_y,
        *(refractive_index_map.get()));
  }

  space_type GetTransmissityTri(
      const Vec3<space_type>& point_on_tri,
      const Tri<space_type,color_type>& tri) {
    if (!transmissity_uv_map_tri) {
      return transmissity_const;
    }

    return ExtractScalarFromTextureMap(
        transmissity_uv_map_tri.value(),
        tri,
        point_on_tri,
        *(transmissity_map.get()));
  }
  space_type GetTransmissitySphere(
      const Vec3<space_type>& point_on_surface,
      const Vec3<space_type>& sphere_origin) {
    if (!transmissity_uv_map_sphere) {
      return transmissity_const;
    }

    return ExtractScalarFromTextureMap(
        point_on_surface,
        sphere_origin,
        transmissity_uv_map_sphere.value().e_y,
        *(transmissity_map.get()));
  }
};

int32_t MapIntoPeriodicBoundary(int32_t x, int32_t size) {
  return x >= 0 ? x % size : size - (x % size);
}

template<typename T, Float_t space_type, Float_t color_type>
T MapTextureNearestNeighbour(
    const Texture<T,space_type,color_type>& map,
    const Vec2<space_type>& uv) {
  int32_t u_discrete = 
      MapIntoPeriodicBoundary(
          static_cast<int32_t>(uv.x * map.im.get()->width),
          map.im.get()->width);
  int32_t v_discrete = 
      MapIntoPeriodicBoundary(
          static_cast<int32_t>(uv.y * map.im.get()->height),
          map.im.get()->height);
  
  return map.im.get()->pixel_ptr[u_discrete + (v_discrete * map.im.get()->width)];
}

//TODO: bilin and bicubuc interpolation also with vec3 and pixel sup

//For Triangles
template<Scalar_t T, Float_t space_type, Float_t color_type>
T ExtractScalarFromTextureMap(
    const UvMapTri<space_type>& UV_coords,
    Tri<space_type,color_type> tri,
    Vec3<space_type> point_on_tri,
    const Texture<T,space_type,color_type>& map) {
  //Roate triangle plane and point so that the normal faces in y direction
  tri.RotateTriToNormalUp();
  point_on_tri.PerformAlignRotation(tri.geometry_normal);
  Vec2<space_type> uv_mapped(point_on_tri.x, point_on_tri.z);

  //UV mapping only allows rotation and uniform scaling (i.e. for both u and v
  //) -> compute these from UV coords and apply transform to point
  Vec2<space_type> edge1_UV = UV_coords.uv_verts[1] - UV_coords.uv_verts[0];
  Vec2<space_type> edge1 =
      Vec2<space_type>(tri.vertex[1].x, tri.vertex[1].z) -
      Vec2<space_type>(tri.vertex[0].x, tri.vertex[0].z);
  Vec2<space_type> offset = UV_coords.uv_verts[0];
  space_type scale = edge1_UV.GetLength() / edge1.GetLength();
  space_type rotation = std::acos(edge1_UV * edge1);

  uv_mapped = (uv_mapped - offset) * scale;
  uv_mapped.RotateCartesian(rotation);
    
  return MapTextureNearestNeighbour(map, uv_mapped);
}
template<Float_t space_type, Float_t color_type>
Pixel<color_type> ExtractPixelFromTextureMap(
    const UvMapTri<space_type>& UV_coords,
    Tri<space_type,color_type> tri,
    Vec3<space_type> point_on_tri,
    const Texture<Pixel<color_type>,space_type,color_type>& map) {
  //Roate triangle plane and point so that the normal faces in y direction
  tri.RotateTriToNormalUp();
  point_on_tri.PerformAlignRotation(tri.geometry_normal);
  Vec2<space_type> uv_mapped(point_on_tri.x, point_on_tri.z);

  //UV mapping only allows rotation and uniform scaling (i.e. for both u and v
  //) -> compute these from UV coords and apply transform to point
  Vec2<space_type> edge1_UV = UV_coords.uv_verts[1] - UV_coords.uv_verts[0];
  Vec2<space_type> edge1 =
      Vec2<space_type>(tri.vertex[1].x, tri.vertex[1].z) -
      Vec2<space_type>(tri.vertex[0].x, tri.vertex[0].z);
  Vec2<space_type> offset = UV_coords.uv_verts[0];
  space_type scale = edge1_UV.GetLength() / edge1.GetLength();
  space_type rotation = std::acos(edge1_UV * edge1);

  uv_mapped = (uv_mapped - offset) * scale;
  uv_mapped.RotateCartesian(rotation);
    
  return MapTextureNearestNeighbour(map, uv_mapped);
}
template<Float_t space_type, Float_t color_type>
Vec3<space_type> ExtractVec3FromTextureMap(
    const UvMapTri<space_type>& UV_coords,
    Tri<space_type,color_type> tri,
    Vec3<space_type> point_on_tri,
    const Texture<Vec3<space_type>,space_type,color_type>& map) {
  //Roate triangle plane and point so that the normal faces in y direction
  tri.RotateTriToNormalUp();
  point_on_tri.PerformAlignRotation(tri.geometry_normal);
  Vec2<space_type> uv_mapped(point_on_tri.x, point_on_tri.z);

  //UV mapping only allows rotation and uniform scaling (i.e. for both u and v
  //) -> compute these from UV coords and apply transform to point
  Vec2<space_type> edge1_UV = UV_coords.uv_verts[1] - UV_coords.uv_verts[0];
  Vec2<space_type> edge1 =
      Vec2<space_type>(tri.vertex[1].x, tri.vertex[1].z) -
      Vec2<space_type>(tri.vertex[0].x, tri.vertex[0].z);
  Vec2<space_type> offset = UV_coords.uv_verts[0];
  space_type scale = edge1_UV.GetLength() / edge1.GetLength();
  space_type rotation = std::acos(edge1_UV * edge1);

  uv_mapped = (uv_mapped - offset) * scale;
  uv_mapped.RotateCartesian(rotation);
    
  return MapTextureNearestNeighbour(map, uv_mapped);
}

//For spheres
template<Scalar_t T, Float_t space_type, Float_t color_type>
T ExtractScalarFromTextureMap(
    Vec3<space_type> point_on_surface,
    const Vec3<space_type>& sphere_origin,
    Vec3<space_type> sphere_pole_axis_normalized,
    const Texture<T,space_type,color_type>& map) {
  //compute phi, theta ; add phi, theta trafo back (r is known)

  //Compute point on "texture-sphere" given the rotation implied by the axis
  sphere_pole_axis_normalized.ConvertToSphericalCoords();
  point_on_surface -= sphere_origin;
  point_on_surface.ConvertToSphericalCoords();
  Vec3<space_type> texture_point =
      (point_on_surface +
       Vec3<space_type>(0.0,sphere_pole_axis_normalized.x,sphere_pole_axis_normalized.z));
  texture_point.ConvertToCartisianCoords();
    
  Vec2<space_type> uv_mapped(0.5 + std::atan2(texture_point.x, texture_point.z)
                          / (2*M_PI),
                    0.5 + std::asin(texture_point.y) / M_PI);
  return MapTextureNearestNeighbour(map, uv_mapped);
}
template<Float_t space_type, Float_t color_type>
Pixel<color_type> ExtractPixelFromTextureMap(
    Vec3<space_type> point_on_surface,
    const Vec3<space_type>& sphere_origin,
    Vec3<space_type> sphere_pole_axis_normalized,
    const Texture<Pixel<color_type>,space_type,color_type>& map) {
  //compute phi, theta ; add phi, theta trafo back (r is known)

  //Compute point on "texture-sphere" given the rotation implied by the axis
  sphere_pole_axis_normalized.ConvertToSphericalCoords();
  point_on_surface -= sphere_origin;
  point_on_surface.ConvertToSphericalCoords();
  Vec3<space_type> texture_point =
      (point_on_surface +
       Vec3<space_type>(0.0,sphere_pole_axis_normalized.x,sphere_pole_axis_normalized.z));
  texture_point.ConvertToCartisianCoords();
    
  Vec2<space_type> uv_mapped(0.5 + std::atan2(texture_point.x, texture_point.z)
                          / (2*M_PI),
                    0.5 + std::asin(texture_point.y) / M_PI);
  return MapTextureNearestNeighbour(map, uv_mapped);
}
template<Float_t space_type, Float_t color_type>
Vec3<space_type> ExtractVec3FromTextureMap(
    Vec3<space_type> point_on_surface,
    const Vec3<space_type>& sphere_origin,
    Vec3<space_type> sphere_pole_axis_normalized,
    const Texture<Vec3<space_type>,space_type,color_type>& map) {
  //compute phi, theta ; add phi, theta trafo back (r is known)

  //Compute point on "texture-sphere" given the rotation implied by the axis
  sphere_pole_axis_normalized.ConvertToSphericalCoords();
  point_on_surface -= sphere_origin;
  point_on_surface.ConvertToSphericalCoords();
  Vec3<space_type> texture_point =
      (point_on_surface +
       Vec3<space_type>(0.0,sphere_pole_axis_normalized.x,sphere_pole_axis_normalized.z));
  texture_point.ConvertToCartisianCoords();
    
  Vec2<space_type> uv_mapped(0.5 + std::atan2(texture_point.x, texture_point.z)
                          / (2*M_PI),
                    0.5 + std::asin(texture_point.y) / M_PI);
  return MapTextureNearestNeighbour(map, uv_mapped);
}