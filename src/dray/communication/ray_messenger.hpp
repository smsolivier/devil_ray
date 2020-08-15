#ifndef DRAY_RAY_MESSENGER_HPP
#define DRAY_RAY_MESSENGER_HPP

#include <mpi.h>
#include <list>
#include <vector>
#include <set>
#include <map>

#include <dray/ray.hpp>
#include <dray/communication/ray_result.hpp>
#include <dray/communication/comm_data.hpp>
#include <dray/communication/messenger.hpp>

// debuggin
#include <iostream>
#include <fstream>
#include <sstream>

#define LOG_MESSAGES 1

namespace dray
{

class MemStream;

class RayMessenger : public Messenger
{
  public:
    RayMessenger(MPI_Comm comm);
    ~RayMessenger();

    void register_messages(int msgSize,
                           int numMsgRecvs,
                           int numICRecvs);

    void send_rays(std::map<int, std::vector<Ray>> &ray_map);
    void send_rays(int dst, std::vector<Ray> &rays);
    bool recv_rays(std::vector<Ray> &rays);

    void send_results(int dst, std::vector<RayResult> &results);
    bool recv_results(std::vector<RayResult> &results);

    // Send/Recv messages.
    void send_msg(int dst, std::vector<int> &msg);
    void send_all_msg(std::vector<int> &msg);
    bool recv_msg(std::vector<MsgCommData> &msgs);

    bool recv_any(std::vector<MsgCommData> *msgs,
                  std::vector<Ray> *rays,
                  std::vector<RayResult> *results,
                  bool blockAndWait);

  private:

#ifdef LOG_MESSAGES
    static int m_message_id;
    std::ofstream m_log;
#endif

    enum
    {
      MESSAGE_TAG = 0xbadbeef,
      RAY_TAG = 0xfeebdab,
      RESULT_TAG = 0xbadbad
    };

    ////Message headers.
    //typedef struct
    //{
    //    int rank, id, tag, numPackets, packet, packetSz, dataSz;
    //} Header;
};

} //namespace dray
#endif
