#include <iostream>
#include <string.h>
#include <dray/communication/memstream.hpp>
#include <dray/communication/ray_messenger.hpp>

using namespace std;

namespace dray
{
#ifdef LOG_MESSAGES
int RayMessenger::m_message_id = 0;
#endif

RayMessenger::RayMessenger(MPI_Comm comm)
  : Messenger(comm)
{
#ifdef LOG_MESSAGES
  std::stringstream ss;
  ss<<"msg_log_"<<this->m_rank<<".txt";
  m_log.open(ss.str());
#endif
}

RayMessenger::~RayMessenger()
{
#ifdef LOG_MESSAGES
  m_log.close();
#endif
}

void
RayMessenger::register_messages(int mSz,
                                int nMsgRecvs,
                                int nICRecvs)
{
    // Msgs are handled as vector<int>.
    // Serialization of msg consists: size_t (num elements) +
    // sender rank + message size.
    int msgSize = sizeof(size_t);
    msgSize += sizeof(int); // sender rank.
    msgSize += (mSz * sizeof(int));

    //During particle advection, the IC state is only serialized.
    int slSize = 256; // avoids splitting send buffer into multiple chunks
    //slSize = sizeof(Ray);
    //slsPerRecv = 64;
    //slsPerRecv = 640; //works
    int slsPerRecv = 10000;

    int dsSize = 2 * sizeof(int);

    this->register_tag(RayMessenger::MESSAGE_TAG, nMsgRecvs, msgSize);
    this->register_tag(RayMessenger::RAY_TAG, nICRecvs, slSize * slsPerRecv);

    this->init_buffers();
}

void
RayMessenger::send_msg(int dst, vector<int> &msg)
{
    MemStream *buff = new MemStream;

    //Write data.
    dray::write(*buff, m_rank);
    dray::write(*buff, msg);
#ifdef LOG_MESSAGES
    dray::write(*buff, m_message_id);
    m_log<<m_rank<<" "<<m_message_id<<"\n";
    m_message_id++;
#endif
    send_data(dst, RayMessenger::MESSAGE_TAG, buff);
//    MsgCnt.value++;
//    CommTime.value += visitTimer->StopTimer(timerHandle, "send_msg");
}

void
RayMessenger::send_all_msg(vector<int> &msg)
{
  for (int i = 0; i < m_mpi_size; i++)
  {
    if (i != m_rank)
    {
      //DBG("          ***************** send_msg to "<<i<<" "<<msg<<endl);
      send_msg(i, msg);
    }
  }
}

bool
RayMessenger::recv_any(vector<MsgCommData> *msgs,
                      std::vector<Ray> *rays,
                      std::vector<RayResult> *results,
                      bool blockAndWait)
{
  set<int> tags;
  if(msgs)
  {
    tags.insert(RayMessenger::MESSAGE_TAG);
    msgs->resize(0);
  }
  if(rays)
  {
    tags.insert(RayMessenger::RAY_TAG);
    rays->resize(0);
  }
  if(results)
  {
    tags.insert(RayMessenger::RESULT_TAG);
    results->resize(0);
  }

  if (tags.empty())
      return false;

  vector<pair<int, MemStream *> > buffers;
  if (! recv_data(tags, buffers, blockAndWait))
      return false;

//    int timerHandle = visitTimer->StartTimer();

  for (size_t i = 0; i < buffers.size(); i++)
  {
      if (buffers[i].first == RayMessenger::MESSAGE_TAG)
      {
        int sendRank;
        vector<int> m;
        dray::read(*buffers[i].second, sendRank);
        dray::read(*buffers[i].second, m);
#ifdef LOG_MESSAGES
        int message_id;
        dray::read(*buffers[i].second, message_id);
        m_log<<sendRank<<" "<<message_id<<"\n";
#endif

        MsgCommData msg(sendRank, m);

        msgs->push_back(msg);
      }
      else if (buffers[i].first == RayMessenger::RAY_TAG)
      {
        int num, sendRank;
        dray::read(*buffers[i].second, sendRank);
        dray::read(*buffers[i].second, num);
#ifdef LOG_MESSAGES
        int message_id;
        dray::read(*buffers[i].second, message_id);
        m_log<<sendRank<<" "<<message_id<<"\n";
#endif
        // TODO: loop through before and allocate mem
        // TODO: read entire buffer at once
        //rays->resize(num + rays->size());
        //DPRINT("["<<m_rank<<"] <-- ["<<sendRank<<"] "<<num<<"\n");
        //std::cout<<"["<<m_rank<<"] <-- ["<<sendRank<<"] ";
        for (int j = 0; j < num; j++)
        {
          Ray ray;
          dray::read(*(buffers[i].second), ray);
          rays->push_back(ray);
          //std::cout<<(*rays)[j].m_pixel_id<<" ";
        }
        //std::cout<<"\n";
      }
      else if (buffers[i].first == RayMessenger::RESULT_TAG)
      {
        int num, sendRank;
        dray::read(*buffers[i].second, sendRank);
        dray::read(*buffers[i].second, num);
#ifdef LOG_MESSAGES
        int message_id;
        dray::read(*buffers[i].second, message_id);
        m_log<<sendRank<<" "<<message_id<<"\n";
#endif
       // TODO: loop through before and allocate mem
       // OR use host memory pool!
       // TODO: read entire buffer at once
       //rays->resize(num + rays->size());
       //DPRINT("["<<m_rank<<"] <-- ["<<sendRank<<"] RES "<<num<<"\n");
       //std::cout<<"["<<m_rank<<"] <-- ["<<sendRank<<"] ";
       for (int j = 0; j < num; j++)
       {
         RayResult result;
         dray::read(*(buffers[i].second), result);
         results->push_back(result);
         //std::cout<<(*rays)[j].m_pixel_id<<" ";
       }
       //std::cout<<"\n";
     }

     delete buffers[i].second;
 }

//    CommTime.value += visitTimer->StopTimer(timerHandle, "recv_any");
    return true;
}

bool
RayMessenger::recv_msg(vector<MsgCommData> &msgs)
{
  return recv_any(&msgs, NULL, NULL, false);
}

void RayMessenger::send_rays(int dst, std::vector<Ray> &rays)
{
  if (dst == m_rank)
  {
      cerr<<"Error. Sending IC to yourself"<<endl;
      return;
  }
  if (rays.empty())
      return;

  MemStream *buff = new MemStream;
  dray::write(*buff, m_rank);
  const int num = rays.size();
  dray::write(*buff, num);

#ifdef LOG_MESSAGES
  dray::write(*buff, m_message_id);
  m_log<<m_rank<<" "<<m_message_id<<"\n";
  m_message_id++;
#endif

    for (auto &ray : rays)
    {
        //std::cout<<"writing "<<ray<<"\n";
        dray::write(*buff, ray);
    }
    send_data(dst, RayMessenger::RAY_TAG, buff);
    rays.clear();
}

void RayMessenger::send_results(int dst, std::vector<RayResult> &results)
{
  if (dst == m_rank)
  {
      cerr<<"Error. Sending result to yourself"<<endl;
      return;
  }
  if (results.empty())
      return;

  MemStream *buff = new MemStream;
  dray::write(*buff, m_rank);
  const int num = results.size();
  dray::write(*buff, num);

#ifdef LOG_MESSAGES
  dray::write(*buff, m_message_id);
  m_log<<m_rank<<" "<<m_message_id<<"\n";
  m_message_id++;
#endif

  for(auto &result: results)
  {
    dray::write(*buff, result);
  }
  send_data(dst, RayMessenger::RESULT_TAG, buff);
  results.clear();
}

void RayMessenger::send_rays(std::map<int, std::vector<Ray>> &ray_map)
{
  for (auto mit = ray_map.begin(); mit != ray_map.end(); mit++)
  {
    if (! mit->second.empty())
    {
      send_rays(mit->first, mit->second);
    }
  }
}

bool RayMessenger::recv_rays(std::vector<Ray> &rays)
{
  return recv_any(NULL, &rays, NULL, false);
}

bool RayMessenger::recv_results(std::vector<RayResult> &results)
{
  return recv_any(NULL, NULL, &results, false);
}

} // namespace dray
