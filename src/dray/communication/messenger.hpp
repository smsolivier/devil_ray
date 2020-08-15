#ifndef DRAY_MESSENGER_HPP
#define DRAY_MESSENGER_HPP

#include <mpi.h>
#include <list>
#include <vector>
#include <set>
#include <map>
//#include "CommData.hpp"

namespace dray
{

class MemStream;

class Messenger
{
  public:
    Messenger(MPI_Comm comm);
    virtual ~Messenger() {}

    //Message headers.
    typedef struct
    {
        int rank, id, tag, numPackets, packet, packetSz, dataSz;
    } Header;

    // Register message tags for this messenger
    // Must be called before init_buffers
    void register_tag(int tag,         // Unique message tag
                      int num_recvs,   // number of receives to check each time
                      int size);       // size in bytes for each message

    // Creates receives buffers for all tags registered to this messenger
    void init_buffers();

    void cleanup()
    {
      cleanup_requests();
    }


    //Manage communication.
    void cleanup_requests(int tag=-1);
    void check_pending_send_requests();

  protected:
    void post_recv(int tag);
    void post_recv(int tag, int sz, int src=-1);
    void send_data(int dst, int tag, MemStream *buff);
    bool recv_data(std::set<int> &tags,
                  std::vector<std::pair<int,MemStream *>> &buffers,
                  bool blockAndWait=false);
    bool recv_data(int tag, std::vector<MemStream *> &buffers,
                   bool blockAndWait=false);

    //void add_header(MemStream *buff);
    //void remove_header(MemStream *input, MemStream *header, MemStream *buff);

    void prepare_for_send(int tag, MemStream *buff, std::vector<unsigned char *> &buffList);
    static bool packet_compare(const unsigned char *a, const unsigned char *b);
    void process_received_buffers(std::vector<unsigned char*> &incomingBuffers,
                                std::vector<std::pair<int, MemStream *>> &buffers);

    // Send/Recv buffer management structures.
    typedef std::pair<MPI_Request, int> RequestTagPair;
    typedef std::pair<int, int> RankIdPair;
    typedef std::map<RequestTagPair, unsigned char *>::iterator bufferIterator;
    typedef std::map<RankIdPair, std::list<unsigned char *>>::iterator packetIterator;

    int m_rank;
    int m_mpi_size;
    MPI_Comm m_mpi_comm;

    std::map<RequestTagPair, unsigned char *> m_send_buffers, m_recv_buffers;
    std::map<RankIdPair, std::list<unsigned char *>> m_recv_packets;

    // Maps MPI_TAG to pair(num buffers, data size).
    std::map<int, std::pair<int, int>> m_message_tag_info;

    //int m_num_msg_recvs;
    //int m_num_sl_recvs;
    //int slSize;
    //int slsPerRecv;
    //int msgSize;
    long m_msg_id; // message counter

};

} // namespace dray
#endif
