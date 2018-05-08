/*
 * BMacRicer.h
 *
 *  Created on: May 4, 2018
 *      Author: Binh
 */

#ifndef INET_LINKLAYER_BMACRICER_BMACRICER_H_
#define INET_LINKLAYER_BMACRICER_BMACRICER_H_
#include "inet/physicallayer/contract/packetlevel/IRadio.h"
#include "inet/linklayer/contract/IMACProtocol.h"
#include "inet/linklayer/common/MACAddress.h"
#include "inet/linklayer/base/MACProtocolBase.h"
#include "inet/linklayer/bmacricer/BMacRicerFrame_m.h"

#include <string>
#include <sstream>
#include <vector>
#include <list>
#include <fstream>

namespace inet {

using namespace physicallayer;

/**
 * @brief Implementation of B-MAC (called also Berkeley MAC, Low Power
 * Listening or LPL).
 *
 * The protocol works as follows: each node is allowed to sleep for
 * slotDuration. After waking up, it first checks the channel for ongoing
 * transmissions.
 * If a transmission is catched (a preamble is received), the node stays awake
 * for at most slotDuration and waits for the actual data packet.
 * If a node wants to send a packet, it first sends preambles for at least
 * slotDuration, thus waking up all nodes in its transmission radius and
 * then sends out the data packet. If a mac-level ack is required, then the
 * receiver sends the ack immediately after receiving the packet (no preambles)
 * and the sender waits for some time more before going back to sleep.
 *
 * B-MAC is designed for low traffic, low power communication in WSN and is one
 * of the most widely used protocols (e.g. it is part of TinyOS).
 * The finite state machine of the protocol is given in the below figure:
 *
 * \image html BMACFSM.png "B-MAC Layer - finite state machine"
 *
 * A paper describing this implementation can be found at:
 * http://www.omnet-workshop.org/2011/uploads/slides/OMNeT_WS2011_S5_C1_Foerster.pdf
 *
 * @class BMacLayer
 * @ingroup macLayer
 * @author Anna Foerster
 *
 */
class INET_API BMacRicerLayer: public MACProtocolBase, public IMACProtocol {
private:
    /** @brief Copy constructor is not allowed.
     */
    BMacRicerLayer(const BMacRicerLayer&);
    /** @brief Assignment operator is not allowed.
     */
    BMacRicerLayer& operator=(const BMacRicerLayer&);

public:
    BMacRicerLayer() :
            macQueue(), nbTxDataPackets(0), nbTxPreambles(0), nbRxDataPackets(
                    0), nbRxPreambles(0), nbMissedAcks(0), nbRecvdAcks(0), nbDroppedDataPackets(
                    0), nbTxAcks(0), macState(INIT), resend_data(NULL), ack_timeout(
            NULL), start_bmac(NULL), wakeup(NULL), send_ack(NULL), cca_timeout(
            NULL), ack_tx_over(NULL), send_preamble(NULL), stop_preambles(
            NULL), data_tx_over(NULL), data_timeout(NULL), lastDataPktSrcAddr(), lastDataPktDestAddr(), txAttempts(
                    0), queueLength(0), animation(
                    false), slotDuration(0), bitrate(0), checkInterval(0), txPower(
                    0), useMacAcks(0), maxTxAttempts(0), stats(false) {
    }
    virtual ~BMacRicerLayer();

    /** @brief Initialization of the module and some variables*/

    virtual void initialize(int) override;
    virtual void refreshDisplay() const override;

    /** @brief Delete all dynamically allocated objects of the module*/
    virtual void finish() override;

    /** @brief Handle messages from lower layer */
    virtual void handleLowerPacket(cPacket *) override;

    /** @brief Handle messages from upper layer */
    virtual void handleUpperPacket(cPacket *) override;

    /** @brief Handle self messages such as timers */
    virtual void handleSelfMessage(cMessage *) override;

    /** @brief Handle control messages from lower layer */
    virtual void receiveSignal(cComponent *source, simsignal_t signalID,
            long value, cObject *details) override;

protected:
    typedef std::list<BMacRicerFrame *> MacQueue;
    /** @brief A queue to store packets from upper layer in case another
     packet is still waiting for transmission.*/
    MacQueue macQueue;

    /** @name Different tracked statistics.*/
    /*@{*/
    long nbTxDataPackets = 0;
    long nbTxPreambles = 0;
    long nbRxDataPackets = 0;
    long nbRxPreambles = 0;
    long nbMissedAcks = 0;
    long nbRecvdAcks = 0;
    long nbDroppedDataPackets = 0;
    long nbTxAcks = 0;
    /*@}*/

    FILE *OutFile;
    FILE *pOutFile;
    FILE *matlabFile;
    FILE *errors;

    int nodeID;

    double deployedTime;

    std::string inputFileForecast = "forecast_2_10.txt"; //forecast profile file name
    std::string inputFileReal = "wind_profile_real_"; //real profile file name

    int numFile = 3;
    char convert[10];

    //ARRSES parameters
    double A = 0;
    double M = 0;
    double F = 0;
    const double Beta = 0.2;
    double Alfa = 0;

    double startShutdown = 0;
    double batteryFailurePeriod = 0;
    int nbWakeUps = 0;
    double startWindow = 0;
    double startTime = 0;
    double endTime = 0;
    double startIdle = 0;
    int stepSize = 370; //for step function
    double forecastInterval = 3600; //prediction length  1h=3600  2h=7200 3h=10800 4h=14400 5h=18000 6h=21600 7h=25200 8h=28800

    double Twi = 10;
    const double recoverParam = 300; //=2.5min
    const double maxTwi = 300; //seconds max interval sleep allowed
    const double minTwi = 1; //seconds min interval sleep allowed
    const double shutdownTwi = 310; //ShutDown interval
    const double beaconTime = 0.05; //period to send a beacon for base station = 50ms
    const double T_idle = 0.052; //period to send a beacon for base station = 52ms
    int K = 10; //wake up times in a slot
    int dimPower_real = 2016; //num o values in real profile
    double P_REAL[2016]; //2016 for Week ---- 8928 for 31days
    int forecastValues = 168; //num o values in forecast profile
    double P_FORECAST[168]; // 1 value for each hour//  168 for Week ----  744 for for 31days

    double E_HARV = 0;
    double e_C = 0; //energy consumed each wake-up interval
    double e_C_n = 0; //energy consumed in a window/slot
    double V_supCap = 4000; //4V initial voltage
    const double Vref = 5000; //Voltage at which to keep constantly supCap
    const double Vmin = 1800; //Voltage at which node must shutdown
    const double Vmax = 5500; //max voltage supported by supercapacitor
    double P_H_r = 0; //real
    double P_H_p = 0; //prevision

    //LUT energy cons PowWow
    const double eta = 0.85; // DC/DC efficiency
    const double C_S = 0.9; //0.9; //0.9F capacitance of supercap
    const double P_LEAK = 43; //43;//microJ
    const double E_WAB = 51; //microJ
    const double E_CCA = 18; //microJ
    const double E_CBT = 9.7; //microJ
    const double E_DT = 80; //80microJ energy cons to transmit pkt
    const double E_ACK = 51; //microJ
    const double P_SLEEP = 85.5; //microW
    const double P_RX = 76.89; //milliW
    const double P_TX = 63.33; //milliW

    /** @brief The MAC address of the interface. */
    MACAddress address;

    /** @brief The radio. */
    IRadio *radio = nullptr;
    IRadio::TransmissionState transmissionState =
            IRadio::TRANSMISSION_STATE_UNDEFINED;

    /** @brief MAC states
     *
     *  The MAC states help to keep track what the MAC is actually
     *  trying to do.
     *  INIT -- node has just started and its status is unclear
     *  SLEEP -- node sleeps, but accepts packets from the network layer
     *  CCA -- Clear Channel Assessment - MAC checks
     *         whether medium is busy
     *  SEND_PREAMBLE -- node sends preambles to wake up all nodes
     *  WAIT_DATA -- node has received at least one preamble from another node
     *                  and wiats for the actual data packet
     *  SEND_DATA -- node has sent enough preambles and sends the actual data
     *                  packet
     *  WAIT_TX_DATA_OVER -- node waits until the data packet sending is ready
     *  WAIT_ACK -- node has sent the data packet and waits for ack from the
     *                 receiving node
     *  SEND_ACK -- node send an ACK back to the sender
     *  WAIT_ACK_TX -- node waits until the transmission of the ack packet is
     *                    over
     */
    enum States {
        INIT,    //0
        SLEEP,    //1
        CCA,    //2
        SEND_PREAMBLE,    //3
        WAIT_DATA,    //4
        SEND_DATA,    //5
        WAIT_TX_DATA_OVER,    //6
        WAIT_ACK,    //7
        SEND_ACK,    //8
        WAIT_ACK_TX,    //9

        RICER_INIT,   //10

        // Rx states
        Rx_WAIT_DATA,   //11
        Rx_BEACON,   //12
        Rx_ACK,   //13
        Rx_BEACON_SENT,   //14
        RICER_Rx_SEND_BEACON,   //15
        RICER_Rx_WAKE_UP,   //

        // Tx states
        Tx_IDLE_LISTENING,   //
        Tx_WUB_RECEIVED,   //
        Tx_CCA,   //
        Tx_CBT,   //
        Tx_FINISH_DATA,   //
        Tx_ACK_RECEIVED,   //
        Tx_SLEEP,   //
        Tx_SHUT_DOWN
    };
    /** @brief The current state of the protocol */
//    States macState = (States) -1;
    States macState;
    /** @brief Types of messages (self messages and packets) the node can
     * process **/
    enum TYPES {
        // packet types
        BMAC_PREAMBLE = 191,
        BMAC_DATA,
        BMAC_ACK,
        // self message types
        BMAC_RESEND_DATA,
        BMAC_ACK_TIMEOUT,
        BMAC_START_BMAC,
        BMAC_WAKE_UP,
        BMAC_SEND_ACK,
        BMAC_CCA_TIMEOUT,
        BMAC_ACK_TX_OVER,
        BMAC_SEND_PREAMBLE,
        BMAC_STOP_PREAMBLES,
        BMAC_DATA_TX_OVER,
        BMAC_DATA_TIMEOUT,

        RICER_BEACON, //205
        RICER_DATA, //206
        RICER_ACK, //207
        RICER_STOP_BEACON, //208
        RICER_STOP_CCA, //209
        RICER_STOP_CBT, //210
        RICER_TIME_OUT, //211
        RICER_SEND_BEACON, //212
        RICER_WAKE_UP,
        RICER_IDLE_TIMEOUT,

        RICER_START,
        RICER_TX_FINISH,
        RICER_WAKEUP,
        RICER_ACK_TX_OVER,
        RICER_BEACON_TX_OVER,
        RICER_SHUTDOWN_WAKEUP
    };

    // messages used in the FSM
    cMessage *resend_data = nullptr;
    cMessage *ack_timeout = nullptr;
    cMessage *start_bmac = nullptr;
    cMessage *wakeup = nullptr;
    cMessage *send_ack = nullptr;
    cMessage *cca_timeout = nullptr;
    cMessage *ack_tx_over = nullptr;
    cMessage *send_preamble = nullptr;
    cMessage *stop_preambles = nullptr;
    cMessage *data_tx_over = nullptr;
    cMessage *data_timeout = nullptr;

    // messages used in the FSM

    cMessage* ricer_start;
    cMessage* ricer_time_out;
    cMessage* ricer_tx_finish;
    cMessage* ricer_wakeup;
    cMessage* ricer_ack;
    cMessage* ricer_ack_tx_over;
    cMessage* ricer_wake_up;

    cMessage* ricer_stop_beacon;
    cMessage* ricer_stop_cca;
    cMessage* ricer_stop_cbt;
    cMessage* ricer_send_beacon;
    cMessage* ricer_idle_timeout;
    cMessage* ricer_shutdown_wakeup;

    cMessage* ricer_beacon_tx_over;

    /** @name Help variables for the acknowledgment process. */
    /*@{*/
    MACAddress lastDataPktSrcAddr;
    MACAddress lastDataPktDestAddr;
    int txAttempts = 0;
    /*@}*/

    /** @brief The maximum length of the queue */
    unsigned int queueLength = 0;
    /** @brief Animate (colorize) the nodes.
     *
     * The color of the node reflects its basic status (not the exact state!)
     * BLACK - node is sleeping
     * GREEN - node is receiving
     * YELLOW - node is sending
     */
    bool animation = false;
    /** @brief The duration of the slot in secs. */
    double slotDuration = 0;
    /** @brief Length of the header*/
    int headerLength = 0;
    /** @brief The bitrate of transmission */
    double bitrate = 0;
    /** @brief Transmission power of the node */
    double txPower;
    /** @brief The duration of CCA */
    double checkInterval = 0;
    /** @brief Use MAC level acks or not */
    bool useMacAcks = false;
    /** @brief Maximum transmission attempts per data packet, when ACKs are
     * used */
    int maxTxAttempts = 0;
    /** @brief Gather stats at the end of the simulation */
    bool stats = false;

    /** @brief Possible colors of the node for animation */
    enum BMAC_COLORS {
        GREEN = 1, BLUE = 2, RED = 3, BLACK = 4, YELLOW = 5
    };

    /** @brief Internal function to change the color of the node */
    void changeDisplayColor(BMAC_COLORS color);

    /** @brief Generate new interface address*/
    virtual void initializeMACAddress();
    virtual InterfaceEntry *createInterfaceEntry() override;
    virtual void handleCommand(cMessage *msg) {
    }

    /** @brief Internal function to send the first packet in the queue */
    void sendDataPacket();

    /** @brief Internal function to send an ACK */
    void sendMacAck();

    /** @brief Internal function to send an BEACON */
    void sendMacBeacon();

    /** @brief Internal function to send one preamble */
    void sendPreamble();

    /** @brief Internal function to attach a signal to the packet */
    void attachSignal(BMacRicerFrame *macPkt);

    /** @brief Internal function to add a new packet from upper to the queue */
    bool addToQueue(cMessage *msg);

    virtual void flushQueue();

    virtual void clearQueue();

    cPacket *decapsMsg(BMacRicerFrame *msg);
    BMacRicerFrame *encapsMsg(cPacket *netwPkt);
    cObject *setUpControlInfo(cMessage * const pMsg,
            const MACAddress& pSrcAddr);

    /************************************** MY FUNCTIONS *********************************/
    /** @brief Internal function to save FORECAST profile */
    void readFileForecast();

    /** @brief Internal function to save real profile */
    void readFileReal();

    /** @brief DEBUG and RESULTS*/
    void writeMatFile();
    void writeErrors(double dt, double P_prev);

    /** @brief Internal function to compute Voltage */
    double calVoltage();

    /** @brief Internal function to compute Voltage when in SHUTDOWN state*/
    double calVoltageShutdown();

    /** @brief Internal function to compute energy harvested */
    double calEnHarv(double tStart, double tEnd);

    /** @brief Internal function to compute energy according to input voltage */
    double calSoC(double voltage);

    /** @brief POWER MANAGER function to compute T_WI */
    void MyPowerManager(double t1, double t2);

    /** @brief Internal function to get forecast value */
    double readForecastValue(double tEnd);

    /** @brief Internal function to compute alfa parameter for ARRSES */
    double calAlfa(double err);

};

} // namespace inet

#endif /* INET_LINKLAYER_BMACRICER_BMACRICER_H_ */
