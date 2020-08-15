#include <iostream>
#include <string.h>
#include <dray/communication/memstream.hpp>
#include <dray/communication/messenger.hpp>
//#include "DebugMeowMeow.hpp"

using namespace std;

namespace dray
{

Messenger::Messenger(MPI_Comm comm)
  : m_mpi_comm(comm)
{
    MPI_Comm_size(comm, &m_mpi_size);
    MPI_Comm_rank(comm, &m_rank);
    m_msg_id = 0;
}

void
Messenger::register_tag(int tag, int num_recvs, int size)
{
  if(m_message_tag_info.find(tag) != m_message_tag_info.end())
  {
    std::cout<<"Warning tag "<<tag<<" already registerd. Overriting\n";
  }

  m_message_tag_info[tag] = pair<int,int>(num_recvs, size);
}

void
Messenger::init_buffers()
{
  //Setup receive buffers.
  map<int, pair<int, int> >::const_iterator it;
  for (it = m_message_tag_info.begin(); it != m_message_tag_info.end(); it++)
  {
    int tag = it->first, num = it->second.first;
    for (int i = 0; i < num; i++)
    {
      post_recv(tag);
    }
  }
}

void
Messenger::cleanup_requests(int tag)
{
  vector<RequestTagPair> delKeys;
  for (bufferIterator i = m_recv_buffers.begin(); i != m_recv_buffers.end(); i++)
  {
    if (tag == -1 || tag == i->first.second)
        delKeys.push_back(i->first);
  }

  if (! delKeys.empty())
  {
    vector<RequestTagPair>::const_iterator it;
    for (it = delKeys.begin(); it != delKeys.end(); it++)
    {
      RequestTagPair v = *it;

      unsigned char *buff = m_recv_buffers[v];
      MPI_Cancel(&(v.first));
      delete [] buff;
      m_recv_buffers.erase(v);
    }
  }
}

void
Messenger::post_recv(int tag)
{
    map<int, pair<int, int> >::const_iterator it = m_message_tag_info.find(tag);
    if (it != m_message_tag_info.end())
    {
      post_recv(tag, it->second.second);
    }
}

void
Messenger::post_recv(int tag, int sz, int src)
{
    sz += sizeof(Messenger::Header);
    unsigned char *buff = new unsigned char[sz];
    memset(buff, 0, sz);

    MPI_Request req;
    if (src == -1)
    {
      MPI_Irecv(buff, sz, MPI_BYTE, MPI_ANY_SOURCE, tag, m_mpi_comm, &req);
    }
    else
    {
      MPI_Irecv(buff, sz, MPI_BYTE, src, tag, m_mpi_comm, &req);
    }

    RequestTagPair entry(req, tag);
    m_recv_buffers[entry] = buff;

    //cerr<<"post_recv: ("<<req<<", "<<tag<<") buff= "<<(void*)buff<<" sz= "<<sz<<endl;
}

void
Messenger::check_pending_send_requests()
{
    bufferIterator it;
    vector<MPI_Request> req, copy;
    vector<int> tags;

    for (it = m_send_buffers.begin(); it != m_send_buffers.end(); it++)
    {
        req.push_back(it->first.first);
        copy.push_back(it->first.first);
        tags.push_back(it->first.second);
    }

    if (req.empty())
        return;

    //See if any sends are done.
    int num = 0, *indices = new int[req.size()];
    MPI_Status *status = new MPI_Status[req.size()];
    int err = MPI_Testsome(req.size(), &req[0], &num, indices, status);
    if (err != MPI_SUCCESS)
    {
        cerr << "Err with MPI_Testsome in PARIC algorithm" << endl;
    }
    for (int i = 0; i < num; i++)
    {
        MPI_Request r = copy[indices[i]];
        int tag = tags[indices[i]];

        RequestTagPair k(r,tag);
        bufferIterator entry = m_send_buffers.find(k);
        if (entry != m_send_buffers.end())
        {
            delete [] entry->second;
            m_send_buffers.erase(entry);
        }
    }

    delete [] indices;
    delete [] status;
}

bool
Messenger::packet_compare(const unsigned char *a, const unsigned char *b)
{
    Messenger::Header ha, hb;
    memcpy(&ha, a, sizeof(ha));
    memcpy(&hb, b, sizeof(hb));

    return ha.packet < hb.packet;
}

void
Messenger::prepare_for_send(int tag,
                            MemStream *buff,
                            vector<unsigned char *> &buffList)
{
  map<int, pair<int, int> >::const_iterator it = m_message_tag_info.find(tag);
  if (it == m_message_tag_info.end())
      throw "message tag not found";

  int bytesLeft = buff->len();
  int maxDataLen = it->second.second;
  Messenger::Header header;
  header.tag = tag;
  header.rank = m_rank;
  header.id = m_msg_id;
  header.numPackets = 1;
  if (buff->len() > (unsigned int)maxDataLen)
      header.numPackets += buff->len() / maxDataLen;

  header.packet = 0;
  header.packetSz = 0;
  header.dataSz = 0;
  m_msg_id++;

  buffList.resize(header.numPackets);
  size_t pos = 0;
  for (int i = 0; i < header.numPackets; i++)
  {
    header.packet = i;
    if (i == (header.numPackets-1))
    {
      header.dataSz = bytesLeft;
    }
    else
    {
      header.dataSz = maxDataLen;
    }

    header.packetSz = header.dataSz + sizeof(header);
    unsigned char *b = new unsigned char[header.packetSz];

    //Write the header.
    unsigned char *bPtr = b;
    memcpy(bPtr, &header, sizeof(header));
    bPtr += sizeof(header);

    //Write the data.
    memcpy(bPtr, &buff->data()[pos], header.dataSz);
    pos += header.dataSz;

    buffList[i] = b;
    bytesLeft -= maxDataLen;
  }

}

void
Messenger::send_data(int dst, int tag, MemStream *buff)
{
  vector<unsigned char *> bufferList;

  //Add headers, break into multiple buffers if needed.
  prepare_for_send(tag, buff, bufferList);

  Messenger::Header header;
  for (size_t i = 0; i < bufferList.size(); i++)
  {
    memcpy(&header, bufferList[i], sizeof(header));
    MPI_Request req;
    int err = MPI_Isend(bufferList[i], header.packetSz, MPI_BYTE, dst,
                        tag, m_mpi_comm, &req);
    if (err != MPI_SUCCESS)
    {
        cerr << "Err with MPI_Isend in PARIC algorithm" << endl;
    }
    //BytesCnt.value += header.packetSz;

    //Add it to m_send_buffers
    RequestTagPair entry(req, tag);
    m_send_buffers[entry] = bufferList[i];
  }

  delete buff;
}

bool
Messenger::recv_data(int tag, std::vector<MemStream *> &buffers,
                            bool blockAndWait)
{
    std::set<int> setTag;
    setTag.insert(tag);
    std::vector<std::pair<int, MemStream *> > b;
    buffers.resize(0);
    if (recv_data(setTag, b, blockAndWait))
    {
        buffers.resize(b.size());
        for (size_t i = 0; i < b.size(); i++)
            buffers[i] = b[i].second;
        return true;
    }
    return false;
}

bool
Messenger::recv_data(set<int> &tags,
                    vector<pair<int, MemStream *> > &buffers,
                    bool blockAndWait)
{
    buffers.resize(0);

    //Find all recv of type tag.
    vector<MPI_Request> req, copy;
    vector<int> reqTags;
    for (bufferIterator i = m_recv_buffers.begin(); i != m_recv_buffers.end(); i++)
    {
        if (tags.find(i->first.second) != tags.end())
        {
            req.push_back(i->first.first);
            copy.push_back(i->first.first);
            reqTags.push_back(i->first.second);
        }
    }

    if (req.empty())
        return false;

    MPI_Status *status = new MPI_Status[req.size()];
    int *indices = new int[req.size()], num = 0;
    if (blockAndWait)
        MPI_Waitsome(req.size(), &req[0], &num, indices, status);
    else
        MPI_Testsome(req.size(), &req[0], &num, indices, status);

    if (num == 0)
    {
        delete [] status;
        delete [] indices;
        return false;
    }

    vector<unsigned char *> incomingBuffers(num);
    for (int i = 0; i < num; i++)
    {
        RequestTagPair entry(copy[indices[i]], reqTags[indices[i]]);
        bufferIterator it = m_recv_buffers.find(entry);
        if ( it == m_recv_buffers.end())
        {
            delete [] status;
            delete [] indices;
            throw "receive buffer not found";
        }

        incomingBuffers[i] = it->second;
        m_recv_buffers.erase(it);
    }

    process_received_buffers(incomingBuffers, buffers);

    for (int i = 0; i < num; i++)
        post_recv(reqTags[indices[i]]);

    delete [] status;
    delete [] indices;

    return ! buffers.empty();
}

void
Messenger::process_received_buffers(vector<unsigned char*> &incomingBuffers,
                                    vector<pair<int, MemStream *> > &buffers)
{
  for (size_t i = 0; i < incomingBuffers.size(); i++)
  {
    unsigned char *buff = incomingBuffers[i];

    //Grab the header.
    Messenger::Header header;
    memcpy(&header, buff, sizeof(header));

    //Only 1 packet, strip off header and add to list.
    if (header.numPackets == 1)
    {
      MemStream *b = new MemStream(header.dataSz, (buff + sizeof(header)));
      b->rewind();
      pair<int, MemStream*> entry(header.tag, b);
      buffers.push_back(entry);
      delete [] buff;
    }
    //Multi packet....
    else
    {
      RankIdPair k(header.rank, header.id);
      packetIterator i2 = m_recv_packets.find(k);

      //First packet. Create a new list and add it.
      if (i2 == m_recv_packets.end())
      {
        list<unsigned char *> l;
        l.push_back(buff);
        m_recv_packets[k] = l;
      }
      else
      {
        i2->second.push_back(buff);

        // The last packet came in, merge into one MemStream.
        if (i2->second.size() == (size_t)header.numPackets)
        {
          //Sort the packets into proper order.
          i2->second.sort(Messenger::packet_compare);

          MemStream *mergedBuff = new MemStream;
          list<unsigned char *>::iterator listIt;

          for (listIt = i2->second.begin(); listIt != i2->second.end(); listIt++)
          {
            unsigned char *bi = *listIt;

            Messenger::Header header;
            memcpy(&header, bi, sizeof(header));
            mergedBuff->write_binary((bi+sizeof(header)), header.dataSz);
            delete [] bi;
          }

          mergedBuff->rewind();
          pair<int, MemStream*> entry(header.tag, mergedBuff);
          buffers.push_back(entry);
          m_recv_packets.erase(i2);
        }
      }
    }
  }
}

}// namespace dray
