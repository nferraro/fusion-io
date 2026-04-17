#ifdef PCMS_ENABLED

#ifndef PCMS_LIB_H
#define PCMS_LIB_H

#include <pcms/coupler.h>
#include <pcms/field.h>
#include <pcms/types.h>
#include <Omega_h_library.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_comm.hpp>

using pcms::Coupler;

namespace {
  template<typename T>
  using remove1Pointer = typename std::remove_pointer<T>::type;
  template<typename T>
  using remove2Pointers = remove1Pointer<remove1Pointer<T>>;
  template<typename T>
  using remove3Pointers = remove1Pointer<remove1Pointer<remove1Pointer<T>>>;
  template <typename T, typename ReturnT>
  using HasOnePointer = typename std::enable_if<std::is_pointer<T>::value 
                                            && !std::is_pointer<remove1Pointer<T>>::value, ReturnT>::type;
  template <typename T, typename ReturnT>
  using HasTwoPointers = typename std::enable_if<std::is_pointer<remove1Pointer<T>>::value 
                                            && !std::is_pointer<remove2Pointers<T>>::value, ReturnT>::type;
  template <typename T, typename ReturnT>
  using HasThreePointers = typename std::enable_if<std::is_pointer<remove2Pointers<T>>::value 
                                              && !std::is_pointer<remove3Pointers<T>>::value, ReturnT>::type;
}

class PCMS_Library {
  public:
  std::unique_ptr<Coupler> coupler;
  pcms::Application* app;

  PCMS_Library(int argc, char** argv, MPI_Comm comm) {
    int isServer = atoi(argv[2]);
    const auto dim = 1;
    auto ranks = redev::LOs({0});
    auto cuts = redev::Reals({0});
    auto ptn = redev::Partition{redev::RCBPtn(dim,ranks,cuts)};
    coupler = std::unique_ptr<Coupler>(new Coupler("pcms_coupler", comm, isServer, ptn));
    app = coupler->AddApplication("pcms_app");
  }
};

template <typename fieldType>
class FusionIOFieldAdapter {
  public:
    using value_type = remove3Pointers<fieldType>;
    using memory_space = pcms::HostMemorySpace;

    FusionIOFieldAdapter(fieldType fieldIn, int sizeX=1, int sizeY=1, int sizeZ=1) :
      field(fieldIn), 
      size{sizeX, sizeY, sizeZ} {
      totalSize = sizeX*sizeY*sizeZ;
      gids = std::vector<pcms::GO>(totalSize);
      std::iota(gids.begin(), gids.end(), static_cast<pcms::GO>(0));
    }

    [[nodiscard]] std::vector<pcms::GO> GetGids() const { 
      return gids;
    }

    [[nodiscard]] pcms::ReversePartitionMap GetReversePartitionMap(const pcms::Partition& partition) const {
      pcms::ReversePartitionMap reverse_partition;
      int local_index = 0;
      for (int i=0; i < totalSize; i++) {
        auto dr = partition.GetDr(i, 0);
        reverse_partition[dr].emplace_back(local_index++);
      }
      return reverse_partition;
    }

    //Serialize 1D Pointers 
    template <typename T = fieldType>
    HasOnePointer<T, int> Serialize(pcms::ScalarArrayView<value_type, memory_space> buffer,
                                    pcms::ScalarArrayView<const pcms::LO, memory_space> permutation) const
    {
      REDEV_ALWAYS_ASSERT(buffer.size() == permutation.size());
      if (buffer.size() > 0) {
        for (int i = 0; i < totalSize; i++) {
          buffer[i] = field[permutation[i]];
        }
      }
      return totalSize;
    }

    //Deserialize 1D Pointers
    template <typename T = fieldType>
    HasOnePointer<T, void> Deserialize(pcms::ScalarArrayView<const value_type, memory_space> buffer,
                                      pcms::ScalarArrayView<const pcms::LO, memory_space> permutation) const
    {
      REDEV_ALWAYS_ASSERT(buffer.size() == permutation.size());
      for (int i = 0; i < buffer.size(); ++i) {
        field[permutation[i]] = buffer[i];
      }
    }

    //Serialize 2D Pointers
    template <typename T = fieldType>
    HasTwoPointers<T, int> Serialize(pcms::ScalarArrayView<value_type, memory_space> buffer,
                  pcms::ScalarArrayView<const pcms::LO, memory_space> permutation) const
    {
      REDEV_ALWAYS_ASSERT(buffer.size() == permutation.size());
      if (buffer.size() > 0) {
        for (int i = 0; i < totalSize; i++) {
          int x = permutation[i] % size[0];
          int y = permutation[i] / size[0];
          buffer[i] = field[x][y];
        }
      }
      return totalSize;
    }

    //Deserialize 2D Pointers
    template <typename T = fieldType>
    HasTwoPointers<T, void> Deserialize(pcms::ScalarArrayView<const value_type, memory_space> buffer,
                                        pcms::ScalarArrayView<const pcms::LO, memory_space> permutation) const
    {
      REDEV_ALWAYS_ASSERT(buffer.size() == permutation.size());
      for (int i = 0; i < buffer.size(); ++i) {
        int x = permutation[i] % size[0];
        int y = permutation[i] / size[0];
        field[x][y] = buffer[i];
      }
    }

    //Serialize 3D Pointers
    template <typename T = fieldType>
    HasThreePointers<T, int> Serialize(pcms::ScalarArrayView<value_type, memory_space> buffer,
                  pcms::ScalarArrayView<const pcms::LO, memory_space> permutation) const
    {
      REDEV_ALWAYS_ASSERT(buffer.size() == permutation.size());
      if (buffer.size() > 0) {
        for (int i = 0; i < totalSize; i++) {
          int z = permutation[i] / (size[0] * size[1]);
          int y = (permutation[i] / size[0]) % size[1];
          int x = permutation[i] % size[0];
          buffer[i] = field[x][y][z];
        }
      }
      return totalSize;
    }

    //Deserialize 3D Pointers
    template <typename T = fieldType>
    HasThreePointers<T, void> Deserialize(pcms::ScalarArrayView<const value_type, memory_space> buffer,
                      pcms::ScalarArrayView<const pcms::LO, memory_space> permutation) const
    {
      REDEV_ALWAYS_ASSERT(buffer.size() == permutation.size());
      for (int i = 0; i < buffer.size(); ++i) {
        int z = permutation[i] / (size[0] * size[1]);
        int y = (permutation[i] / size[0]) % size[1];
        int x = permutation[i] % size[0];
        field[x][y][z] = buffer[i];
      }
    }

  private:
    fieldType field;
    int size[3];
    int totalSize;
    std::vector<pcms::GO> gids;

}; //class FusionIOFieldAdapter

#endif //PCMS_LIB_H

#endif //PCMS_ENABLED