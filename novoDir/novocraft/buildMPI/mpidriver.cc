//
// File:   mpi_driver.cc
// Author: Colin Hercus
// Copyright 2009 Colin Hercus
//
// Created on 26 Nov 2009
//

#include <stdio.h>
#include "mpidriver.h"
#include <mpi.h>

const unsigned char mpi_driver::eof_message[4] = "EOF";
const NV_Datatype mpi_driver::eof_datatype = DT_CHAR;
MPI_Datatype dtConvert[DT_SIZE] = {MPI_CHAR, MPI_INT, MPI_UNSIGNED};

struct msgbuffers {
    int size;
    MPI_Request *requests;
    char **buffer;
    unsigned int *buffer_size;
    msgbuffers()
    {
        requests = NULL;
        buffer = NULL;
        buffer_size = NULL;
        size = 0;
    }

    ~msgbuffers()
    {
        mpitracel;
        delete[] requests;
        for(int i = 0; i < size; i++) {
            mpitracef("%d", i);
            if(buffer[i] != NULL)
                delete[] buffer[i];
        }
        delete[] buffer_size;
        delete[] buffer;
        mpitracef("Exit");
    }

    void allocate(int n, unsigned int bufsize)
    {
        mpitracef("%d buffers of %d bytes", n, bufsize);
        size = n;
        requests = new MPI_Request[size];
        buffer_size = new unsigned int[size];
        buffer = new char *[size];
        for(int i = 0; i < size; i++)
        {
            requests[i] = MPI_REQUEST_NULL;
            buffer_size[i] = bufsize;
            buffer[i] = new char[buffer_size[i]];
        }
    }

    void resize(int i, unsigned int bufsize) {
        if(bufsize < buffer_size[i]) return;
        mpitracef("%d %d", i, bufsize);
        if(buffer[i] != NULL)
            delete[] buffer[i];
        buffer_size[i] = bufsize;
        buffer[i] = new char[buffer_size[i]];
    }
};


mpi_driver::mpi_driver() {};

mpi_driver::~mpi_driver()
{
    mpitracef("destruct %d", rank);
}

//Starts MPI Process and decides whether we are a master or slave
void mpi_driver::run_mpi(int argc, char** argv)
{
    int rc;
    int got_threading;
    rc = MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &got_threading);
    mpitracef("mpi_driver::run_mpi() Threading %d got %d", MPI_THREAD_FUNNELED, got_threading);
    if (rc != MPI_SUCCESS) {
        fprintf (stderr, "Error initialising %s, MPI::Init failed. Terminating.\n", argv[0]);
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    mpitracel;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    mpitracel;
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    if (num_tasks <= 1) {
        fprintf (stderr, "Error initialising %s, no slave tasks. Terminating.\n", argv[0]);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    mpitracef("Rank %d, Num %d", rank, num_tasks);
    if(rank != 0)
        run_slaves(argc, argv);
    else
        run_master(argc, argv);

    MPI_Finalize();

}

void mpi_driver::run_master(int argc, char** argv)
{
    static const int SENDBUFFERSPERSLAVE = 2;
    static const int RECVBUFFERS = num_tasks;
    int task, buf_idx;
    unsigned int tag;
    int block_num = 0;
    long numMessages = 0;
    masterRequest = new msgbuffers;

    mpitracel;
    init_master(argc, argv);
    mpitracel;
    masterRequest->allocate(RECVBUFFERS + num_tasks*SENDBUFFERSPERSLAVE, max_length());
    eof_sent = new bool[num_tasks];
    wait_sync = new bool[num_tasks];
    for(int i = 0; i < num_tasks; i++)
        eof_sent[i] = false;
    mpitracel;
    int num_running;
    for(int i = 1; i < RECVBUFFERS; i++)
        MPI_Irecv(masterRequest->buffer[i], masterRequest->buffer_size[i], MPI_CHAR, i, MPI_ANY_TAG, MPI_COMM_WORLD, &(masterRequest->requests)[i]);
    while(!eof()) {
        num_running = num_tasks - 1;
        for(int i = 0; i < num_tasks; i++)
            wait_sync[i] = false;
        mpitracel;
        for(buf_idx = RECVBUFFERS; buf_idx < masterRequest->size ; buf_idx++) {
            task = (buf_idx - RECVBUFFERS) % num_tasks;
            if(task == 0) continue;
            mpitracef("task %d buffer %d", task, buf_idx);
            int msg_sz = next_message(task, masterRequest->buffer[buf_idx], masterRequest->buffer_size[buf_idx], tag);
            MPI_Issend(masterRequest->buffer[buf_idx], msg_sz, MPI_CHAR, task, tag, MPI_COMM_WORLD, &masterRequest->requests[buf_idx]);
        }
        mpitracel;

        while(num_running > 0)
        {
            int flag, req_idx;
            MPI_Status status;
            mpitracef("Testany");
            MPI_Waitany(masterRequest->size, masterRequest->requests, &buf_idx, &status);
            mpitracef("testany %d %d %d", buf_idx, num_tasks, num_running);
            if(buf_idx == MPI_UNDEFINED)
                continue;
            if(buf_idx >= RECVBUFFERS)
            {
                task = (buf_idx - RECVBUFFERS) % num_tasks;
                numMessages++;
                if(!eof() && !syncPoint())
                {
                    mpitracef("send message %d", task);
                    int msg_sz = next_message(task, masterRequest->buffer[buf_idx], masterRequest->buffer_size[buf_idx], tag);
                    mpitracef("Issend task %d buf %d len %d", task, buf_idx, msg_sz);
                    MPI_Issend(masterRequest->buffer[buf_idx], msg_sz, MPI_CHAR, task, tag, MPI_COMM_WORLD, &masterRequest->requests[buf_idx]);
                    continue;
                }
                if(eof()) {
                    if(!eof_sent[task]) {
                        mpitracef("send eof %d", task);
                        MPI_Issend((void *)eof_message, eof_len, dtConvert[eof_datatype], task, TAG_EOF, MPI_COMM_WORLD, &masterRequest->requests[buf_idx]);
                        eof_sent[task] = true;
                    }
                }
                else if(syncPoint() && !wait_sync[task]) {
                    mpitracef("send SYNC %d", task);
                    MPI_Issend((void *)eof_message, 0, dtConvert[eof_datatype], task, TAG_SYNC, MPI_COMM_WORLD, &masterRequest->requests[buf_idx]);
                    wait_sync[task] = true;
                }
                continue;
            }
            else
            {
                numMessages++;
                mpitracef("Master Recvd %d from %d", status.MPI_TAG, status.MPI_SOURCE);
                int i = buf_idx;
                int slave = status.MPI_SOURCE;
                if(status.MPI_TAG == TAG_MSGSIZE) {
                    unsigned int newsize = *(unsigned int *)masterRequest->buffer;
                    mpitracef("From:%d resize buf[%d] %d", slave, i, newsize);
                    masterRequest->resize(i, newsize);
                    MPI_Irecv(masterRequest->buffer[i], masterRequest->buffer_size[i], MPI_CHAR, i, MPI_ANY_TAG, MPI_COMM_WORLD, &masterRequest->requests[i]);
                    continue;
                }
                if(status.MPI_TAG == TAG_SYNCD) {
                    mpitracef("SYNCD from %d Running %d Buf[%d]", slave, num_running, i);
                    MPI_Irecv(masterRequest->buffer[i], masterRequest->buffer_size[i], MPI_CHAR, i, MPI_ANY_TAG, MPI_COMM_WORLD, &masterRequest->requests[i]);
                    num_running--;
                    continue;
                }
                if(status.MPI_TAG == TAG_EOF) {
                    mpitracef("EOF from %d Running %d Buf[%d]", slave, num_running, i);
                    if(num_running > RECVBUFFERS)
                        MPI_Irecv(masterRequest->buffer[i], masterRequest->buffer_size[i], MPI_CHAR, i, MPI_ANY_TAG, MPI_COMM_WORLD, &masterRequest->requests[i]);
                    num_running--;
                    continue;
                }
                int msg_len;
                (void)MPI_Get_count(&status, MPI_CHAR, &msg_len);
                master_process_message(slave, masterRequest->buffer[buf_idx], status.MPI_TAG, msg_len);
                MPI_Irecv(masterRequest->buffer[i], masterRequest->buffer_size[i], MPI_CHAR, i, MPI_ANY_TAG, MPI_COMM_WORLD, &masterRequest->requests[i]);
                continue;
            }
        }
        if(!eof())
            synced();
    }
    mpitracef("Barrier");
    MPI_Barrier(MPI_COMM_WORLD);
    final_results();
    mpitracef("exit");

}

void mpi_driver::init_master(int argc, char** argv) {};

bool mpi_driver::final_results() { return true; };

void mpi_driver::run_slaves(int argc, char** argv)
{
    mpitracel;
    init_slave(argc, argv);
    mpitracel;
    eof_received = false;
//    MPI_Request request = MPI_REQUEST_NULL;
    MPI_Status status;
    while(true)
    {
        unsigned char eof_buf[eof_len];
        mpitracef("R%d wait Probe()", rank);
        MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        mpitracef("R%d probe success", rank);
        int tag = status.MPI_TAG;
        if(tag == TAG_EOF) {
            mpitracef("R%d eof recvd", rank);
            MPI_Recv(eof_buf, mpi_driver::eof_len,  dtConvert[mpi_driver::eof_datatype], 0, TAG_EOF, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            break;
        }
        mpitracef("R%d",rank);
        int msg_len;
        (void)MPI_Get_count(&status, MPI_CHAR, &msg_len);
        process_message(tag, msg_len);
    }
    eof_slave();
    mpitracef("Barrier - Slave %d", rank);
    MPI_Barrier(MPI_COMM_WORLD);
    final_slave();
    mpitracef("R%d exit slave", rank);
}

void mpi_driver::init_slave(int argc, char** argv) {};

void mpi_driver::eof_slave() {};

void mpi_driver::final_slave() {};

int NV_MPI_Send(void* buf, int count, NV_Datatype datatype, int task, int tag)
{
    mpitracef("%lx %d %d %d", buf, count, datatype, dtConvert[datatype]);
    return MPI_Send(buf, count, dtConvert[datatype], task, tag, MPI_COMM_WORLD);
}

int NV_MPI_Recv(void * msg, int msg_len,  NV_Datatype datatype)
{
    MPI_Recv(msg, msg_len,  dtConvert[datatype], MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

int NV_MPI_MsgLen(int count, NV_Datatype datatype)
{
    MPI_Aint size;
    MPI_Type_extent(dtConvert[datatype], &size);
    return size * count;
}

void Die()
{
    MPI_Abort(MPI_COMM_WORLD, -1);
}

int NV_MPI_Bcast(void* buf, int count, NV_Datatype datatype)
{
    return MPI_Bcast(buf, count, dtConvert[datatype], 0, MPI_COMM_WORLD);
}

int Barrier()
{
    return MPI_Barrier(MPI_COMM_WORLD);
}
