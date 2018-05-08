/*
 * BMacRicer.cc
 *
 *  Created on: May 4, 2018
 *      Author: Binh
 */

#include "inet/common/INETUtils.h"
#include "inet/common/INETMath.h"
#include "inet/networklayer/common/InterfaceEntry.h"
#include "inet/common/ModuleAccess.h"
#include "inet/linklayer/contract/IMACProtocolControlInfo.h"
#include "inet/linklayer/common/SimpleLinkLayerControlInfo.h"
#include "inet/linklayer/bmacricer/BMacRicer.h"

namespace inet {

Define_Module(BMacRicerLayer);

void BMacRicerLayer::initialize(int stage) {
    EV_DETAIL
                     << "--------------------------------------------------------------------"
                     << endl;
    EV_DETAIL
                     << "--------------------------xoxoxo------------------------------------"
                     << endl;
    MACProtocolBase::initialize(stage);
    if (stage == INITSTAGE_LOCAL) {
        queueLength = par("queueLength");
        queueLength = hasPar("queueLength") ? par("queueLength") : 10;
        animation = hasPar("animation") ? par("animation") : true;
        slotDuration = hasPar("slotDuration") ? par("slotDuration") : 1.;
        bitrate = hasPar("bitrate") ? par("bitrate") : 15360.;
        headerLength = hasPar("headerLength") ? par("headerLength") : 10.;
        checkInterval = hasPar("checkInterval") ? par("checkInterval") : 0.1;
        txPower = hasPar("txPower") ? par("txPower") : 50.;
        useMacAcks = hasPar("useMACAcks") ? par("useMACAcks") : false;
        maxTxAttempts = hasPar("maxTxAttempts") ? par("maxTxAttempts") : 2;

        EV_DETAIL << "headerLength: " << headerLength << ", bitrate: "
                         << bitrate << endl;
        stats = par("stats");
        nbTxDataPackets = 0;
        nbTxPreambles = 0;
        nbRxDataPackets = 0;
        nbRxPreambles = 0;
        nbMissedAcks = 0;
        nbRecvdAcks = 0;
        nbDroppedDataPackets = 0;
        nbTxAcks = 0;

        txAttempts = 0;
        lastDataPktDestAddr = MACAddress::BROADCAST_ADDRESS;
        lastDataPktSrcAddr = MACAddress::BROADCAST_ADDRESS;

        macState = RICER_INIT;

        initializeMACAddress();
        registerInterface();

        cModule *radioModule = getModuleFromPar<cModule>(par("radioModule"),
                this);
        radioModule->subscribe(IRadio::radioModeChangedSignal, this);
        radioModule->subscribe(IRadio::transmissionStateChangedSignal, this);
        radio = check_and_cast<IRadio *>(radioModule);

        // init the dropped packet info
        WATCH(macState);
    } else if (stage == INITSTAGE_LINK_LAYER) {
        wakeup = new cMessage("wakeup");
        wakeup->setKind(BMAC_WAKE_UP);

        data_timeout = new cMessage("data_timeout");
        data_timeout->setKind(BMAC_DATA_TIMEOUT);
        data_timeout->setSchedulingPriority(100);

        data_tx_over = new cMessage("data_tx_over");
        data_tx_over->setKind(BMAC_DATA_TX_OVER);

        stop_preambles = new cMessage("stop_preambles");
        stop_preambles->setKind(BMAC_STOP_PREAMBLES);

        send_preamble = new cMessage("send_preamble");
        send_preamble->setKind(BMAC_SEND_PREAMBLE);

        ack_tx_over = new cMessage("ack_tx_over");
        ack_tx_over->setKind(BMAC_ACK_TX_OVER);

        cca_timeout = new cMessage("cca_timeout");
        cca_timeout->setKind(BMAC_CCA_TIMEOUT);
        cca_timeout->setSchedulingPriority(100);

        send_ack = new cMessage("send_ack");
        send_ack->setKind(BMAC_SEND_ACK);

        start_bmac = new cMessage("start_bmac");
        start_bmac->setKind(BMAC_START_BMAC);

        ack_timeout = new cMessage("ack_timeout");
        ack_timeout->setKind(BMAC_ACK_TIMEOUT);

        resend_data = new cMessage("resend_data");
        resend_data->setKind(BMAC_RESEND_DATA);
        resend_data->setSchedulingPriority(100);

        ricer_start = new cMessage("ricer_start");
        ricer_start->setKind(RICER_START);

        ricer_tx_finish = new cMessage("ricer_tx_finish");
        ricer_tx_finish->setKind(RICER_TX_FINISH);

        ricer_wakeup = new cMessage("ricer_wakeup");
        ricer_wakeup->setKind(RICER_WAKEUP);

        ricer_ack = new cMessage("ricer_ack");
        ricer_ack->setKind(RICER_ACK);

        ricer_ack_tx_over = new cMessage("ricer_ack_tx_over");
        ricer_ack_tx_over->setKind(RICER_ACK_TX_OVER);

        ricer_stop_beacon = new cMessage("ricer_stop_beacon");
        ricer_stop_beacon->setKind(RICER_STOP_BEACON);

        ricer_stop_cca = new cMessage("ricer_stop_cca");
        ricer_stop_cca->setKind(RICER_STOP_CCA);

        ricer_stop_cbt = new cMessage("ricer_stop_cbt");
        ricer_stop_cbt->setKind(RICER_STOP_CBT);

        ricer_beacon_tx_over = new cMessage("ricer_beacon_tx_over");
        ricer_beacon_tx_over->setKind(RICER_BEACON_TX_OVER);

        ricer_time_out = new cMessage("ricer_time_out");
        ricer_time_out->setKind(RICER_TIME_OUT);

        ricer_send_beacon = new cMessage("ricer_send_beacon");
        ricer_send_beacon->setKind(RICER_SEND_BEACON);

        ricer_wake_up = new cMessage("ricer_wake_up");
        ricer_wake_up->setKind(RICER_WAKE_UP);

        ricer_idle_timeout = new cMessage("ricer_idle_timeout");
        ricer_idle_timeout->setKind(RICER_IDLE_TIMEOUT);

        ricer_shutdown_wakeup = new cMessage("ricer_shutdown_wakeup");
        ricer_shutdown_wakeup->setKind(RICER_SHUTDOWN_WAKEUP);
        deployedTime =
                static_cast<double>(this->getAncestorPar("deployedTime"));
        scheduleAt(simTime() + deployedTime, ricer_start);
    }
}

BMacRicerLayer::~BMacRicerLayer() {
    cancelAndDelete(wakeup);
    cancelAndDelete(data_timeout);
    cancelAndDelete(data_tx_over);
    cancelAndDelete(stop_preambles);
    cancelAndDelete(send_preamble);
    cancelAndDelete(ack_tx_over);
    cancelAndDelete(cca_timeout);
    cancelAndDelete(send_ack);
    cancelAndDelete(start_bmac);
    cancelAndDelete(ack_timeout);
    cancelAndDelete(resend_data);

    for (auto & elem : macQueue) {
        delete (elem);
    }
    macQueue.clear();
}

void BMacRicerLayer::finish() {
    recordScalar("nbTxDataPackets", nbTxDataPackets);
    recordScalar("nbTxPreambles", nbTxPreambles);
    recordScalar("nbRxDataPackets", nbRxDataPackets);
    recordScalar("nbRxPreambles", nbRxPreambles);
    recordScalar("nbMissedAcks", nbMissedAcks);
    recordScalar("nbRecvdAcks", nbRecvdAcks);
    recordScalar("nbTxAcks", nbTxAcks);
    recordScalar("nbDroppedDataPackets", nbDroppedDataPackets);
    //recordScalar("timeSleep", timeSleep);
    //recordScalar("timeRX", timeRX);
    //recordScalar("timeTX", timeTX);
}

void BMacRicerLayer::initializeMACAddress() {
    const char *addrstr = par("address");

    if (!strcmp(addrstr, "auto")) {
        // assign automatic address
        address = MACAddress::generateAutoAddress();

        // change module parameter from "auto" to concrete address
        par("address").setStringValue(address.str().c_str());
    } else {
        address.setAddress(addrstr);
    }
}

InterfaceEntry *BMacRicerLayer::createInterfaceEntry() {
    InterfaceEntry *e = new InterfaceEntry(this);

    // data rate
    e->setDatarate(bitrate);

    // generate a link-layer address to be used as interface token for IPv6
    e->setMACAddress(address);
    e->setInterfaceToken(address.formInterfaceIdentifier());

    // capabilities
    e->setMtu(par("mtu").longValue());
    e->setMulticast(false);
    e->setBroadcast(true);

    return e;
}

/**
 * Check whether the queue is not full: if yes, print a warning and drop the
 * packet. Then initiate sending of the packet, if the node is sleeping. Do
 * nothing, if node is working.
 */
void BMacRicerLayer::handleUpperPacket(cPacket *msg) {
    bool pktAdded = addToQueue(msg);
    if (!pktAdded)
        return;
    // force wakeup now
    if (wakeup->isScheduled() && (macState == SLEEP)) {
        cancelEvent(wakeup);
        scheduleAt(simTime() + dblrand() * 0.1f, wakeup);
    }
}

/**
 * Send one short preamble packet immediately.
 */
void BMacRicerLayer::sendPreamble() {
    BMacRicerFrame *preamble = new BMacRicerFrame();
    preamble->setSrcAddr(address);
    preamble->setDestAddr(MACAddress::BROADCAST_ADDRESS);
    preamble->setKind(BMAC_PREAMBLE);
    preamble->setBitLength(headerLength);

    //attach signal and send down
    attachSignal(preamble);
    sendDown(preamble);
    nbTxPreambles++;
}

/**
 * Send one short preamble packet immediately.
 */
void BMacRicerLayer::sendMacAck() {
    BMacRicerFrame *ack = new BMacRicerFrame();
    ack->setSrcAddr(address);
    ack->setDestAddr(lastDataPktSrcAddr);
    ack->setKind(RICER_ACK);
    ack->setBitLength(headerLength);
    ack->setSrcProcId(nodeID);

    //attach signal and send down
    attachSignal(ack);
    sendDown(ack);
    nbTxAcks++;
    //endSimulation();
}

/**
 * Send one BEACON packet immediately.
 */
void BMacRicerLayer::sendMacBeacon() {
    BMacRicerFrame *wub = new BMacRicerFrame();
    wub->setSrcAddr(address);
    wub->setDestAddr(MACAddress::BROADCAST_ADDRESS);
    wub->setKind(RICER_BEACON);
    wub->setBitLength(headerLength);
    wub->setSrcProcId(nodeID);

    attachSignal(wub);
    sendDown(wub);
}

/**
 * Handle BMAC preambles and received data packets.
 */
void BMacRicerLayer::handleLowerPacket(cPacket *msg) {
    if (msg->hasBitError()) {
        EV << "Received " << msg
                  << " contains bit errors or collision, dropping it\n";
        delete msg;
        return;
    } else
        // simply pass the massage as self message, to be processed by the FSM.
        handleSelfMessage(msg);
}

void BMacRicerLayer::sendDataPacket() {
//    nbTxDataPackets++;
//    BMacRicerFrame *pkt = macQueue.front()->dup();
//    attachSignal(pkt);
//    lastDataPktDestAddr = pkt->getDestAddr();
//    pkt->setKind(BMAC_DATA);
//    sendDown(pkt);
    BMacRicerFrame *data = new BMacRicerFrame();
    data->setSrcAddr(address);
    data->setDestAddr(MACAddress::BROADCAST_ADDRESS);
    data->setKind(RICER_DATA);
    data->setBitLength(headerLength);
    data->setSrcProcId(nodeID);

    attachSignal(data);
    sendDown(data);
}

//read file with forecast values
void BMacRicerLayer::readFileForecast() {
    double log_local_wake_up;

    std::ifstream inFile(inputFileForecast);
    std::string line = "";
    int i = 0;
    while (std::getline(inFile, line)) {
        std::istringstream ff(line);
        ff >> log_local_wake_up;
        P_FORECAST[i++] = log_local_wake_up;

        if (P_FORECAST[i - 1] < 140) //140 microW as min threshold
            P_FORECAST[i - 1] = 0;
        if (P_FORECAST[i - 1] > 851) //851 microW as max threshold
            P_FORECAST[i - 1] = 851;
    }
    inFile.close();
}

//calculate voltage in shutDown state(only P_Leak dissipation)
double BMacRicerLayer::calVoltageShutdown() {
    double E_init = calSoC(V_supCap); // 0.5 * C_S * V_supCap * V_supCap; //energy when node deployed(n-1)
    E_init -= (P_LEAK * (endTime - startTime));

    //must add harvested energy
    double E_fin = E_init + calEnHarv(startTime, endTime);

    if (E_fin <= 0)
        E_fin = 0;
    return sqrt(E_fin * 2 / C_S); //end  Voltage
}

//calculate voltage
double BMacRicerLayer::calVoltage() {
    double E_init = calSoC(V_supCap); // 0.5 * C_S * V_supCap * V_supCap; //energy when node deployed(n-1)
    //remove energy consumption
    E_init -= (e_C / eta);
    E_init -= (P_SLEEP * (endTime - startTime) / eta);
    E_init -= (P_LEAK * (endTime - startTime));
    //must add harvested energy
    double E_fin = E_init + calEnHarv(startTime, endTime);

    return sqrt(E_fin * 2 / C_S); //end  Voltage
}

//get forecast power
double BMacRicerLayer::readForecastValue(double tEnd) {
    int i = (int) tEnd / forecastInterval;
    if (i < forecastValues) {
        return P_FORECAST[i];
    }
    return -1;
}

//read file with 5 granularity samples(profile)
void BMacRicerLayer::readFileReal() {
    double log_local_wake_up;
    itoa(numFile, convert, 10);
    inputFileReal.append(convert);
    inputFileReal.append(".txt");

    std::ifstream inFile(inputFileReal);
    std::string line = "";
    int i = 0;
    while (std::getline(inFile, line)) {
        std::istringstream ff(line);
        ff >> log_local_wake_up;
        P_REAL[i++] = log_local_wake_up;

        if (P_REAL[i - 1] < 140) //140 microW as min threshold
            P_REAL[i - 1] = 0;
        if (P_REAL[i - 1] > 851) //851 microW as max threshold
            P_REAL[i - 1] = 851;
    }
    inFile.close();
}

//get harvested energy in a slot n (using profile not forecast!!!)
double BMacRicerLayer::calEnHarv(double tStart, double tEnd) {
    int i = (int) tStart / (300);
    int j = (int) tEnd / (300);
    if (i < dimPower_real && j < dimPower_real) {
        double eH = 0;
        i++;
        j++;

        if (i == j)
            eH = (tEnd - tStart) * P_REAL[i - 1];
        else {
            if (j - i > 1)
                for (int c = 1, d = i; c < j - i; c++, d++) {
                    eH += ((300) * P_REAL[d]);
                }

            eH += (P_REAL[i - 1] * ((300 * i) - tStart));
            eH += (P_REAL[j - 1] * (tEnd - 300 * (j - 1)));
        }
        return eH;
    }
    return 0;
}

double BMacRicerLayer::calSoC(double voltage) {
    if (voltage <= 0)
        return 0;

    return 0.5 * C_S * voltage * voltage;
}

void BMacRicerLayer::MyPowerManager(double t1, double t2) {
    //file output
    fprintf(pOutFile, "\n\nPower Manager is activated, %f\n", endTime);
    double num = 0;
    double den = 0;
    double budgetEnergy = 0;

    double dt = (calEnHarv(t1, t2) / (t2 - t1));
    double P_prev = readForecastValue(t2);
    double P_prev_2h = readForecastValue(t2 + forecastInterval);
    if (P_prev == -1)
        P_prev = dt;
    if (P_prev_2h == -1)
        P_prev_2h = dt;

    Alfa = calAlfa(dt - P_prev);
    double gamma = (t2 / 3600 - (int) (t2 / 3600));
    //PM Policy
#ifdef PAST
    P_H_p = dt;
#endif
#ifdef ARRSES_1
    P_H_p = Alfa*dt + (1-Alfa)*P_prev;
#endif
#ifdef ARRSES_2
    P_H_p = Alfa*dt + (1-Alfa)*((1-gamma)*P_prev + (gamma)*P_prev_2h);
#endif

    budgetEnergy = 0.5 * C_S * (V_supCap * V_supCap - Vref * Vref);

    //battery fully charged
    if (V_supCap >= Vmax) {
        Twi = minTwi;
        V_supCap = Vmax; //cannot have more than Vmax
    } else {
        //calculate numerator and denominator
        num = (e_C_n - budgetEnergy * eta) / K; //P_C_n - budgetEnergy * n ;
        den = eta * (P_H_p - P_LEAK) - P_SLEEP;

        if (den <= 0) {
#ifdef STEP_FN
            Twi = ceil(Twi * ((Vmax - V_supCap)/stepSize));
#else
            Twi = maxTwi;
#endif
        } else {
            double T_WU_n = ceil(num / den);
            Twi = T_WU_n;
        }
        //filter for condition 1 <= Twi <= 300
        if (Twi > maxTwi)
            Twi = maxTwi;
        if (Twi < minTwi)
            Twi = 1;
    }
    fprintf(pOutFile, "**2**Next wake up : %f\n", Twi);
    fprintf(pOutFile, "**3**Actual Voltage : %f\n", V_supCap);
    fprintf(pOutFile, "**4**Budget energy : %f\n", budgetEnergy);
    fprintf(pOutFile, "**5**Total consumed energy : %f\n", e_C_n);
    fprintf(pOutFile, "**6**Real harvested energy : %f\n", E_HARV);
    fprintf(pOutFile, "**7**Prevision harvested power : %f\n", P_H_p);
    fprintf(pOutFile, "**8**tstart tend : %f, %f\n", startTime, endTime);
    fprintf(pOutFile,
            "**9 alfa: %f,prev real PH: %f, forcast 1: %f,prev forcast 2: %f, forcast: %f\n",
            Alfa, dt, P_prev, P_prev_2h, P_H_p);
    fprintf(pOutFile, "**10**num/den : %f\n", ceil(num / den));
    writeErrors(dt, P_prev);

}

//ARRSES algorithm
double BMacRicerLayer::calAlfa(double err) {
    A = Beta * err + (1 - Beta) * A;
    M = Beta * fabs(err) + (1 - Beta) * M;
    if (M == 0)
        M = 1;
    fprintf(pOutFile, "**10 A: %f,M: %f, err: %f, abs: %f\n", A, M, err,
            fabs(A / M));

    return fabs(A / M);
}
void BMacRicerLayer::writeMatFile() {
    fprintf(matlabFile, "%f %f %f %f %f %f %f\n", endTime, Twi, V_supCap, P_H_r,
            E_HARV, batteryFailurePeriod, P_H_p);
}
void BMacRicerLayer::writeErrors(double dt, double P_prev) {
    fprintf(errors, "%f %f %f %f\n", endTime, dt, P_prev, P_H_p);
}

/**
 * Handle own messages:
 * BMAC_WAKEUP: wake up the node, check the channel for some time.
 * BMAC_CHECK_CHANNEL: if the channel is free, check whether there is something
 * in the queue and switch the radio to TX. When switched to TX, the node will
 * start sending preambles for a full slot duration. If the channel is busy,
 * stay awake to receive message. Schedule a timeout to handle false alarms.
 * BMAC_SEND_PREAMBLES: sending of preambles over. Next time the data packet
 * will be send out (single one).
 * BMAC_TIMEOUT_DATA: timeout the node after a false busy channel alarm. Go
 * back to sleep.
 */
void BMacRicerLayer::handleSelfMessage(cMessage *msg) {
    EV << "----------------------------------------------------------------"
              << endl;
    EV << "--------------------------xoxoxo------------------------------------"
              << endl;
    switch (macState) {
    case RICER_INIT:
        if (msg->getKind() == RICER_START) {
            if (nodeID == 0) {
                EV << "****Rx: I am the base station" << endl;
                changeDisplayColor(GREEN);
                radio->setRadioMode(IRadio::RADIO_MODE_TRANSMITTER);
                scheduleAt(simTime(), ricer_send_beacon);
                macState = RICER_Rx_SEND_BEACON;
            } else {
                //open file for results
                pOutFile = fopen("result0.txt", "w");
                matlabFile = fopen("matlabResult.txt", "w");
                errors = fopen("matlabErrors.txt", "w");

                //read forecast and 5 granularity samples
                readFileForecast();
                readFileReal();

                startWindow = startIdle = startTime = simTime().dbl();
                EV << "****Tx: I am the end device..waiting for beacon... "
                          << endl;
                changeDisplayColor(GREEN);
                radio->setRadioMode(IRadio::RADIO_MODE_RECEIVER); //wait for WAB
                scheduleAt(simTime() + T_idle, ricer_idle_timeout); //sleep period
                macState = Tx_IDLE_LISTENING;
            }
            return;
        }
        break;
    case Tx_IDLE_LISTENING:
        if (msg->getKind() == RICER_BEACON) {
            EV << "****Tx: detected beacon.. " << endl;
            cancelEvent(ricer_idle_timeout);
            double idleTime = simTime().dbl() - startIdle;
            e_C += idleTime * P_RX * 1000; // *1000 because PRX is in mW -> idle energy consumption
            macState = Tx_WUB_RECEIVED;
            scheduleAt(simTime() + Twi, ricer_wakeup); //sleep period
            scheduleAt(simTime(), ricer_stop_beacon);
            return;
        }
        if (msg->getKind() == RICER_IDLE_TIMEOUT) {
            changeDisplayColor(BLACK);
            radio->setRadioMode(IRadio::RADIO_MODE_SLEEP);
            scheduleAt(simTime() + Twi, ricer_wakeup); //sleep period
            macState = Tx_SLEEP;
            return;
        }
        break;

    case Tx_WUB_RECEIVED:
        if (msg->getKind() == RICER_STOP_BEACON) {
            EV << "****Tx: Beacon received... " << endl;
            e_C += E_WAB; //WAB energy consumption
            macState = Tx_CCA;
            scheduleAt(simTime(), ricer_stop_cca);
            return;
        }
        break;

    case Tx_CCA:
        if (msg->getKind() == RICER_STOP_CCA) {
            EV << "****Tx: CCA done.. " << endl;
            e_C += E_CCA; //CCA energy consumption
            macState = Tx_CBT;
            scheduleAt(simTime(), ricer_stop_cbt);
            return;
        }
        break;

    case Tx_CBT:
        if (msg->getKind() == RICER_STOP_CBT) {
            EV << "****Tx: CBT done.. " << endl;
            e_C += E_CBT; //CBT idle energy consumption
            macState = Tx_FINISH_DATA;
            //changeDisplayColor(GREEN);
            radio->setRadioMode(IRadio::RADIO_MODE_TRANSMITTER);
            sendDataPacket();
            EV << "****Tx: Sending packet.. " << endl;
            return;
        }
        break;

    case Tx_FINISH_DATA:
        if (msg->getKind() == RICER_TX_FINISH) {
            e_C += E_DT; //DT idle energy consumption
            EV << "****Tx: Data packet is sent..waiting for ACK" << endl;
            changeDisplayColor(YELLOW);
            radio->setRadioMode(IRadio::RADIO_MODE_RECEIVER);
            macState = Tx_ACK_RECEIVED;
            return;
        }
        break;
    case Tx_ACK_RECEIVED:
        if (msg->getKind() == RICER_ACK) {
            e_C += E_ACK; //ACK idle energy consumption
            EV << "****Tx: ACK packet from node n." << msg->getSrcProcId()
                      << endl;
            EV << "****Tx: Going to sleep..." << endl;
            changeDisplayColor(BLACK);
            radio->setRadioMode(IRadio::RADIO_MODE_SLEEP);
            macState = Tx_SLEEP;
            return;
        }

        break;
    case Tx_SLEEP:
        if (msg->getKind() == RICER_WAKEUP) {
            /*check SoC*/
            endTime = simTime().dbl();
            V_supCap = calVoltage();
            E_HARV = calEnHarv(startTime, endTime);
            P_H_r = E_HARV / (endTime - startTime);
            if (V_supCap <= Vmin) {
                startShutdown = endTime;
                Twi = shutdownTwi;
                //reset parameters
                nbWakeUps = 0;
                e_C_n = 0;
                e_C = 0;
                //writeMatFile();
                scheduleAt(simTime() + Twi, ricer_shutdown_wakeup); //sleep period
                macState = Tx_SHUT_DOWN;
            } else {
                e_C_n += e_C;
                nbWakeUps++;
                if (nbWakeUps == K) // 10 times wake up
                        {
                    E_HARV = calEnHarv(startWindow, endTime);
                    P_H_r = E_HARV / (endTime - startWindow);
                    MyPowerManager(startWindow, endTime);
                    EV << "****Tx: Wake up: I compute Energy harvested on time "
                              << startTime << " " << endTime
                              << " equal to new TWI: " << Twi << endl; // P_HARV[endTime/(Twi + deployedTime)] << endl;

                    //reset parameters
                    startWindow = endTime; // startTime = endTime;
                    nbWakeUps = 0;
                    e_C_n = 0;

                    writeMatFile();
                }
                //reset parameter
                e_C = 0;

                startIdle = startTime = endTime;
                changeDisplayColor(GREEN);
                EV << "****Tx: Idle for beacon..." << endl;
                radio->setRadioMode(IRadio::RADIO_MODE_RECEIVER); //wait for WAB
                macState = Tx_IDLE_LISTENING;
            }
            return;
        } else { //simulation is over!!!
            fclose(pOutFile);
            fclose(matlabFile);
            fclose(errors);
            return;
        }

        break;

    case Tx_SHUT_DOWN:
        if (msg->getKind() == RICER_SHUTDOWN_WAKEUP) {
            endTime = simTime().dbl();
            batteryFailurePeriod += (endTime - startShutdown);
            V_supCap = calVoltageShutdown(); //check out actual voltage if it has increased

            E_HARV = calEnHarv(startTime, endTime);
            P_H_r = E_HARV / (endTime - startTime); //just for tracking

            //I should stay here
            if (V_supCap < Vmin) {
                if (V_supCap <= 0)
                    V_supCap = 0;

                startShutdown = simTime().dbl();
                scheduleAt(simTime() + Twi, ricer_shutdown_wakeup); //sleep period
                macState = Tx_SHUT_DOWN;
            } else {
                Twi = recoverParam; //default after waking up from shutdown recoverParam
                startWindow = startIdle = simTime().dbl();
                changeDisplayColor(GREEN);
                EV << "****Tx: Idle for beacon..." << endl;
                radio->setRadioMode(IRadio::RADIO_MODE_RECEIVER); //wait for WAB
                macState = Tx_IDLE_LISTENING;
            }
            writeMatFile();

            //reset parameters
            startTime = endTime;

            return;
        } else { //simulation is over!!!
            fclose(pOutFile);
            fclose(matlabFile);
            fclose(errors);
            return;
        }

        break;
    case RICER_Rx_SEND_BEACON:
        if (msg->getKind() == RICER_SEND_BEACON) {
            //set the next wake up event
            EV << "****Rx: Sending BEACON.... : " << endl;
            scheduleAt(simTime() + beaconTime, ricer_wake_up); //cur_wake_up/scale_down, ricer_wake_up);
            sendMacBeacon();
            macState = Rx_BEACON_SENT;
            return;
        }
        break;
    case Rx_BEACON_SENT:
        if (msg->getKind() == RICER_BEACON_TX_OVER) //from physical layer
                {
            EV << "****Rx: BEACON packet sent!" << endl;
            changeDisplayColor(YELLOW);
            radio->setRadioMode(IRadio::RADIO_MODE_RECEIVER);
            macState = Rx_WAIT_DATA;
            scheduleAt(simTime() + 0.02, ricer_time_out); //set data timeout to 20ms
            return;
        }
        break;

    case Rx_WAIT_DATA:
        if (msg->getKind() == RICER_TIME_OUT) {
            EV << "****Rx: No data coming..." << endl;
            radio->setRadioMode(IRadio::RADIO_MODE_SLEEP); //turn off radio
            changeDisplayColor(BLACK);
            macState = RICER_Rx_WAKE_UP;

            return;
        }
        if (msg->getKind() == RICER_DATA) {
            EV << "****Rx: I have a data packet from node n."
                      << msg->getSrcProcId() << endl;
            cancelEvent(ricer_time_out);

            EV << "****Rx: Now I answer with an ACK" << endl;
            EV << "****Rx: I am sending a packet..." << endl;
            changeDisplayColor(GREEN);
            radio->setRadioMode(IRadio::RADIO_MODE_TRANSMITTER);
            sendMacAck();
            macState = Rx_ACK;
            return;
        }
        break;
    case Rx_ACK:
        if (msg->getKind() == RICER_ACK_TX_OVER) {
            EV << "****Rx: ACK packet is sent" << endl;
            changeDisplayColor(BLACK);
            EV << "****Rx: ...going to SLEEP" << endl;
            radio->setRadioMode(IRadio::RADIO_MODE_SLEEP);
            macState = RICER_Rx_WAKE_UP;
            return;
        }
        break;

    case RICER_Rx_WAKE_UP:
        if (msg->getKind() == RICER_WAKE_UP) {
            EV << "****Rx: Wake up.." << endl;

            radio->setRadioMode(IRadio::RADIO_MODE_TRANSMITTER);
            changeDisplayColor(YELLOW);
            macState = RICER_Rx_SEND_BEACON;
            scheduleAt(simTime(), ricer_send_beacon);
            return;
        }

        break;
//    case INIT:
//        if (msg->getKind() == BMAC_START_BMAC) {
//            EV_DETAIL << "State INIT, message BMAC_START, new state SLEEP"
//                             << endl;
//            radio->setRadioMode(IRadio::RADIO_MODE_SLEEP);
//            macState = SLEEP;
//            scheduleAt(simTime() + dblrand() * slotDuration, wakeup);
//            return;
//        }
//        break;
//
//    case SLEEP:
//        if (msg->getKind() == BMAC_WAKE_UP) {
//            EV_DETAIL << "State SLEEP, message BMAC_WAKEUP, new state CCA"
//                             << endl;
//            scheduleAt(simTime() + checkInterval, cca_timeout);
//            radio->setRadioMode(IRadio::RADIO_MODE_RECEIVER);
//            macState = CCA;
//            return;
//        }
//        break;
//
//    case CCA:
//        if (msg->getKind() == BMAC_CCA_TIMEOUT) {
//            // channel is clear
//            // something waiting in eth queue?
//            if (macQueue.size() > 0) {
//                EV_DETAIL << "State CCA, message CCA_TIMEOUT, new state"
//                        " SEND_PREAMBLE" << endl;
//                macState = SEND_PREAMBLE;
//                radio->setRadioMode(IRadio::RADIO_MODE_TRANSMITTER);
//                scheduleAt(simTime() + slotDuration, stop_preambles);
//                return;
//            }
//            // if not, go back to sleep and wake up after a full period
//            else {
//                EV_DETAIL << "State CCA, message CCA_TIMEOUT, new state SLEEP"
//                                 << endl;
//                scheduleAt(simTime() + slotDuration, wakeup);
//                macState = SLEEP;
//                radio->setRadioMode(IRadio::RADIO_MODE_SLEEP);
//                return;
//            }
//        }
//        // during CCA, we received a preamble. Go to state WAIT_DATA and
//        // schedule the timeout.
//        if (msg->getKind() == BMAC_PREAMBLE) {
//            nbRxPreambles++;
//            EV_DETAIL << "State CCA, message BMAC_PREAMBLE received, new state"
//                    " WAIT_DATA" << endl;
//            macState = WAIT_DATA;
//            cancelEvent(cca_timeout);
//            scheduleAt(simTime() + slotDuration + checkInterval, data_timeout);
//            delete msg;
//            return;
//        }
//        // this case is very, very, very improbable, but let's do it.
//        // if in CCA and the node receives directly the data packet, switch to
//        // state WAIT_DATA and re-send the message
//        if (msg->getKind() == BMAC_DATA) {
//            nbRxDataPackets++;
//            EV_DETAIL << "State CCA, message BMAC_DATA, new state WAIT_DATA"
//                             << endl;
//            macState = WAIT_DATA;
//            cancelEvent(cca_timeout);
//            scheduleAt(simTime() + slotDuration + checkInterval, data_timeout);
//            scheduleAt(simTime(), msg);
//            return;
//        }
//        //in case we get an ACK, we simply dicard it, because it means the end
//        //of another communication
//        if (msg->getKind() == BMAC_ACK) {
//            EV_DETAIL << "State CCA, message BMAC_ACK, new state CCA" << endl;
//            delete msg;
//            return;
//        }
//        break;
//
//    case SEND_PREAMBLE:
//        if (msg->getKind() == BMAC_SEND_PREAMBLE) {
//            EV_DETAIL << "State SEND_PREAMBLE, message BMAC_SEND_PREAMBLE, new"
//                    " state SEND_PREAMBLE" << endl;
//            sendPreamble();
//            scheduleAt(simTime() + 0.5f * checkInterval, send_preamble);
//            macState = SEND_PREAMBLE;
//            return;
//        }
//        // simply change the state to SEND_DATA
//        if (msg->getKind() == BMAC_STOP_PREAMBLES) {
//            EV_DETAIL << "State SEND_PREAMBLE, message BMAC_STOP_PREAMBLES, new"
//                    " state SEND_DATA" << endl;
//            macState = SEND_DATA;
//            txAttempts = 1;
//            return;
//        }
//        break;
//
//    case SEND_DATA:
//        if ((msg->getKind() == BMAC_SEND_PREAMBLE)
//                || (msg->getKind() == BMAC_RESEND_DATA)) {
//            EV_DETAIL << "State SEND_DATA, message BMAC_SEND_PREAMBLE or"
//                    " BMAC_RESEND_DATA, new state WAIT_TX_DATA_OVER" << endl;
//            // send the data packet
//            sendDataPacket();
//            macState = WAIT_TX_DATA_OVER;
//            return;
//        }
//        break;
//
//    case WAIT_TX_DATA_OVER:
//        if (msg->getKind() == BMAC_DATA_TX_OVER) {
//            if ((useMacAcks) && !lastDataPktDestAddr.isBroadcast()) {
//                EV_DETAIL
//                                 << "State WAIT_TX_DATA_OVER, message BMAC_DATA_TX_OVER,"
//                                         " new state WAIT_ACK" << endl;
//                macState = WAIT_ACK;
//                radio->setRadioMode(IRadio::RADIO_MODE_RECEIVER);
//                scheduleAt(simTime() + checkInterval, ack_timeout);
//            } else {
//                EV_DETAIL
//                                 << "State WAIT_TX_DATA_OVER, message BMAC_DATA_TX_OVER,"
//                                         " new state  SLEEP" << endl;
//                delete macQueue.front();
//                macQueue.pop_front();
//                // if something in the queue, wakeup soon.
//                if (macQueue.size() > 0)
//                    scheduleAt(simTime() + dblrand() * checkInterval, wakeup);
//                else
//                    scheduleAt(simTime() + slotDuration, wakeup);
//                macState = SLEEP;
//                radio->setRadioMode(IRadio::RADIO_MODE_SLEEP);
//            }
//            return;
//        }
//        break;
//
//    case WAIT_ACK:
//        if (msg->getKind() == BMAC_ACK_TIMEOUT) {
//            // No ACK received. try again or drop.
//            if (txAttempts < maxTxAttempts) {
//                EV_DETAIL
//                                 << "State WAIT_ACK, message BMAC_ACK_TIMEOUT, new state"
//                                         " SEND_DATA" << endl;
//                txAttempts++;
//                macState = SEND_PREAMBLE;
//                scheduleAt(simTime() + slotDuration, stop_preambles);
//                radio->setRadioMode(IRadio::RADIO_MODE_TRANSMITTER);
//            } else {
//                EV_DETAIL
//                                 << "State WAIT_ACK, message BMAC_ACK_TIMEOUT, new state"
//                                         " SLEEP" << endl;
//                //drop the packet
//                cMessage *mac = macQueue.front();
//                macQueue.pop_front();
//                emit(NF_LINK_BREAK, mac);
//                delete mac;
//
//                // if something in the queue, wakeup soon.
//                if (macQueue.size() > 0)
//                    scheduleAt(simTime() + dblrand() * checkInterval, wakeup);
//                else
//                    scheduleAt(simTime() + slotDuration, wakeup);
//                macState = SLEEP;
//                radio->setRadioMode(IRadio::RADIO_MODE_SLEEP);
//                nbMissedAcks++;
//            }
//            return;
//        }
//        //ignore and other packets
//        if ((msg->getKind() == BMAC_DATA)
//                || (msg->getKind() == BMAC_PREAMBLE)) {
//            EV_DETAIL
//                             << "State WAIT_ACK, message BMAC_DATA or BMAC_PREMABLE, new"
//                                     " state WAIT_ACK" << endl;
//            delete msg;
//            return;
//        }
//        if (msg->getKind() == BMAC_ACK) {
//            EV_DETAIL << "State WAIT_ACK, message BMAC_ACK" << endl;
//            BMacRicerFrame *mac = static_cast<BMacRicerFrame *>(msg);
//            const MACAddress src = mac->getSrcAddr();
//            // the right ACK is received..
//            EV_DETAIL << "We are waiting for ACK from : " << lastDataPktDestAddr
//                             << ", and ACK came from : " << src << endl;
//            if (src == lastDataPktDestAddr) {
//                EV_DETAIL << "New state SLEEP" << endl;
//                nbRecvdAcks++;
//                lastDataPktDestAddr = MACAddress::BROADCAST_ADDRESS;
//                cancelEvent(ack_timeout);
//                delete macQueue.front();
//                macQueue.pop_front();
//                // if something in the queue, wakeup soon.
//                if (macQueue.size() > 0)
//                    scheduleAt(simTime() + dblrand() * checkInterval, wakeup);
//                else
//                    scheduleAt(simTime() + slotDuration, wakeup);
//                macState = SLEEP;
//                radio->setRadioMode(IRadio::RADIO_MODE_SLEEP);
//                lastDataPktDestAddr = MACAddress::BROADCAST_ADDRESS;
//            }
//            delete msg;
//            return;
//        }
//        break;
//
//    case WAIT_DATA:
//        if (msg->getKind() == BMAC_PREAMBLE) {
//            //nothing happens
//            EV_DETAIL << "State WAIT_DATA, message BMAC_PREAMBLE, new state"
//                    " WAIT_DATA" << endl;
//            nbRxPreambles++;
//            delete msg;
//            return;
//        }
//        if (msg->getKind() == BMAC_ACK) {
//            //nothing happens
//            EV_DETAIL
//                             << "State WAIT_DATA, message BMAC_ACK, new state WAIT_DATA"
//                             << endl;
//            delete msg;
//            return;
//        }
//        if (msg->getKind() == BMAC_DATA) {
//            nbRxDataPackets++;
//            BMacRicerFrame *mac = static_cast<BMacRicerFrame *>(msg);
//            const MACAddress& dest = mac->getDestAddr();
//            const MACAddress& src = mac->getSrcAddr();
//            if ((dest == address) || dest.isBroadcast()) {
//                EV_DETAIL << "Local delivery " << mac << endl;
//                sendUp(decapsMsg(mac));
//            } else {
//                EV_DETAIL << "Received " << mac
//                                 << " is not for us, dropping frame." << endl;
//                delete msg;
//                msg = nullptr;
//                mac = nullptr;
//            }
//
//            cancelEvent(data_timeout);
//            if ((useMacAcks) && (dest == address)) {
//                EV_DETAIL << "State WAIT_DATA, message BMAC_DATA, new state"
//                        " SEND_ACK" << endl;
//                macState = SEND_ACK;
//                lastDataPktSrcAddr = src;
//                radio->setRadioMode(IRadio::RADIO_MODE_TRANSMITTER);
//            } else {
//                EV_DETAIL
//                                 << "State WAIT_DATA, message BMAC_DATA, new state SLEEP"
//                                 << endl;
//                // if something in the queue, wakeup soon.
//                if (macQueue.size() > 0)
//                    scheduleAt(simTime() + dblrand() * checkInterval, wakeup);
//                else
//                    scheduleAt(simTime() + slotDuration, wakeup);
//                macState = SLEEP;
//                radio->setRadioMode(IRadio::RADIO_MODE_SLEEP);
//            }
//            return;
//        }
//        if (msg->getKind() == BMAC_DATA_TIMEOUT) {
//            EV_DETAIL << "State WAIT_DATA, message BMAC_DATA_TIMEOUT, new state"
//                    " SLEEP" << endl;
//            // if something in the queue, wakeup soon.
//            if (macQueue.size() > 0)
//                scheduleAt(simTime() + dblrand() * checkInterval, wakeup);
//            else
//                scheduleAt(simTime() + slotDuration, wakeup);
//            macState = SLEEP;
//            radio->setRadioMode(IRadio::RADIO_MODE_SLEEP);
//            return;
//        }
//        break;
//
//    case SEND_ACK:
//        if (msg->getKind() == BMAC_SEND_ACK) {
//            EV_DETAIL << "State SEND_ACK, message BMAC_SEND_ACK, new state"
//                    " WAIT_ACK_TX" << endl;
//            // send now the ack packet
//            sendMacAck();
//            macState = WAIT_ACK_TX;
//            return;
//        }
//        break;
//
//    case WAIT_ACK_TX:
//        if (msg->getKind() == BMAC_ACK_TX_OVER) {
//            EV_DETAIL
//                             << "State WAIT_ACK_TX, message BMAC_ACK_TX_OVER, new state"
//                                     " SLEEP" << endl;
//            // ack sent, go to sleep now.
//            // if something in the queue, wakeup soon.
//            if (macQueue.size() > 0)
//                scheduleAt(simTime() + dblrand() * checkInterval, wakeup);
//            else
//                scheduleAt(simTime() + slotDuration, wakeup);
//            macState = SLEEP;
//            radio->setRadioMode(IRadio::RADIO_MODE_SLEEP);
//            lastDataPktSrcAddr = MACAddress::BROADCAST_ADDRESS;
//            return;
//        }
//        break;
    }
    throw cRuntimeError(
            "Undefined event of type %d in state %d (radio mode %d, radio reception state %d, radio transmission state %d)!",
            msg->getKind(), macState, radio->getRadioMode(),
            radio->getReceptionState(), radio->getTransmissionState());
}

void BMacRicerLayer::receiveSignal(cComponent *source, simsignal_t signalID,
        long value, cObject *details) {
    EV_DETAIL
                     << "------------------------123-----------------------------------------"
                     << endl;
    EV_DETAIL
                     << "--------------------------xoxoxo------------------------------------"
                     << endl;
    Enter_Method_Silent();
    if (signalID == IRadio::radioModeChangedSignal) {
//        IRadio::RadioMode radioMode = (IRadio::RadioMode)value;
//        if (radioMode == IRadio::RADIO_MODE_TRANSMITTER) {
//            // we just switched to TX after CCA, so simply send the first
//            // sendPremable self message
//            if (macState == SEND_PREAMBLE)
//            scheduleAt(simTime(), send_preamble);
//            else if (macState == SEND_ACK)
//            scheduleAt(simTime(), send_ack);
//            // we were waiting for acks, but none came. we switched to TX and now
//            // need to resend data
//            else if (macState == SEND_DATA)
//            scheduleAt(simTime(), resend_data);
//
//        }
        IRadio::RadioMode radioMode = (IRadio::RadioMode) value;
        if ((macState == SEND_PREAMBLE)
                && (radioMode == IRadio::RADIO_MODE_TRANSMITTER)) {
            scheduleAt(simTime(), send_preamble);
        }
        if ((macState == SEND_ACK)
                && (radioMode == IRadio::RADIO_MODE_TRANSMITTER)) {
            scheduleAt(simTime(), send_ack);
        }
        // we were waiting for acks, but none came. we switched to TX and now
        // need to resend data
        if ((macState == SEND_DATA)
                && (radioMode == IRadio::RADIO_MODE_TRANSMITTER)) {
            scheduleAt(simTime(), resend_data);
        }
    }
    // Transmission of one packet is over
    else if (signalID == IRadio::transmissionStateChangedSignal) {
//        IRadio::TransmissionState newRadioTransmissionState = (IRadio::TransmissionState)value;
//        if (transmissionState == IRadio::TRANSMISSION_STATE_TRANSMITTING && newRadioTransmissionState == IRadio::TRANSMISSION_STATE_IDLE) {
//            if (macState == WAIT_TX_DATA_OVER)
//            scheduleAt(simTime(), data_tx_over);
//            else if (macState == WAIT_ACK_TX)
//            scheduleAt(simTime(), ack_tx_over);
//        }
//        transmissionState = newRadioTransmissionState;
        if (macState == Tx_FINISH_DATA) {
            scheduleAt(simTime(), ricer_tx_finish);
        }
        if (macState == Rx_ACK) {
            scheduleAt(simTime(), ricer_ack_tx_over);
        }
        if (macState == Rx_BEACON_SENT) {
            scheduleAt(simTime(), ricer_beacon_tx_over);
        }
    }
}

/**
 * Encapsulates the received network-layer packet into a BMacFrame and set all
 * needed header fields.
 */
bool BMacRicerLayer::addToQueue(cMessage *msg) {
    if (macQueue.size() >= queueLength) {
        // queue is full, message has to be deleted
        EV_DETAIL << "New packet arrived, but queue is FULL, so new packet is"
                " deleted\n";
        emit(packetFromUpperDroppedSignal, msg);
        nbDroppedDataPackets++;
        return false;
    }

    BMacRicerFrame *macPkt = encapsMsg((cPacket *) msg);
    macQueue.push_back(macPkt);
    EV_DETAIL << "Max queue length: " << queueLength << ", packet put in queue"
            "\n  queue size: " << macQueue.size() << " macState: " << macState
                     << endl;
    return true;
}

void BMacRicerLayer::flushQueue() {
    // TODO:
    macQueue.clear();
}

void BMacRicerLayer::clearQueue() {
    macQueue.clear();
}

void BMacRicerLayer::attachSignal(BMacRicerFrame *macPkt) {
    //calc signal duration
    simtime_t duration = macPkt->getBitLength() / bitrate;
    //create and initialize control info with new signal
    macPkt->setDuration(duration);
}

/**
 * Change the color of the node for animation purposes.
 */
void BMacRicerLayer::refreshDisplay() const {
    if (!animation)
        return;
    cDisplayString& dispStr = findContainingNode(this)->getDisplayString();

    switch (macState) {
    case INIT:
    case SLEEP:
        dispStr.setTagArg("b", 3, "black");
        break;

    case CCA:
        dispStr.setTagArg("b", 3, "green");
        break;

    case SEND_ACK:
    case SEND_PREAMBLE:
    case SEND_DATA:
        dispStr.setTagArg("b", 3, "blue");
        break;

    case WAIT_ACK:
    case WAIT_DATA:
    case WAIT_TX_DATA_OVER:
    case WAIT_ACK_TX:
        dispStr.setTagArg("b", 3, "yellow");
        break;

    default:
        dispStr.setTagArg("b", 3, "");
        break;
    }
}

/*void BMacLayer::changeMacState(States newState)
 {
 switch (macState)
 {
 case RX:
 timeRX += (simTime() - lastTime);
 break;
 case TX:
 timeTX += (simTime() - lastTime);
 break;
 case SLEEP:
 timeSleep += (simTime() - lastTime);
 break;
 case CCA:
 timeRX += (simTime() - lastTime);
 }
 lastTime = simTime();

 switch (newState)
 {
 case CCA:
 changeDisplayColor(GREEN);
 break;
 case TX:
 changeDisplayColor(BLUE);
 break;
 case SLEEP:
 changeDisplayColor(BLACK);
 break;
 case RX:
 changeDisplayColor(YELLOW);
 break;
 }

 macState = newState;
 }*/

cPacket *BMacRicerLayer::decapsMsg(BMacRicerFrame *msg) {
    cPacket *m = msg->decapsulate();
    setUpControlInfo(m, msg->getSrcAddr());
    // delete the macPkt
    delete msg;
    EV_DETAIL << " message decapsulated " << endl;
    return m;
}

BMacRicerFrame *BMacRicerLayer::encapsMsg(cPacket *netwPkt) {
    BMacRicerFrame *pkt = new BMacRicerFrame(netwPkt->getName(),
            netwPkt->getKind());
    pkt->setBitLength(headerLength);

    // copy dest address from the Control Info attached to the network
    // message by the network layer
    IMACProtocolControlInfo *cInfo = check_and_cast<IMACProtocolControlInfo *>(
            netwPkt->removeControlInfo());
    EV_DETAIL << "CInfo removed, mac addr=" << cInfo->getDestinationAddress()
                     << endl;
    pkt->setDestAddr(cInfo->getDestinationAddress());

    //delete the control info
    delete cInfo;

    //set the src address to own mac address (nic module getId())
    pkt->setSrcAddr(address);

    //encapsulate the network packet
    pkt->encapsulate(netwPkt);
    EV_DETAIL << "pkt encapsulated\n";

    return pkt;
}

/**
 * Attaches a "control info" (MacToNetw) structure (object) to the message pMsg.
 */
cObject *BMacRicerLayer::setUpControlInfo(cMessage * const pMsg,
        const MACAddress& pSrcAddr) {
    SimpleLinkLayerControlInfo * const cCtrlInfo =
            new SimpleLinkLayerControlInfo();
    cCtrlInfo->setSrc(pSrcAddr);
    cCtrlInfo->setInterfaceId(interfaceEntry->getInterfaceId());
    pMsg->setControlInfo(cCtrlInfo);
    return cCtrlInfo;
}

void BMacRicerLayer::changeDisplayColor(BMAC_COLORS color) {
    if (!animation)
        return;
    cDisplayString& dispStr = this->getDisplayString();
    //b=40,40,rect,black,black,2"
    if (color == GREEN)
        dispStr.setTagArg("b", 3, "green");
    //dispStr.parse("b=40,40,rect,green,green,2");
    if (color == BLUE)
        dispStr.setTagArg("b", 3, "blue");
    //dispStr.parse("b=40,40,rect,blue,blue,2");
    if (color == RED)
        dispStr.setTagArg("b", 3, "red");
    //dispStr.parse("b=40,40,rect,red,red,2");
    if (color == BLACK)
        dispStr.setTagArg("b", 3, "black");
    //dispStr.parse("b=40,40,rect,black,black,2");
    if (color == YELLOW)
        dispStr.setTagArg("b", 3, "yellow");
    //dispStr.parse("b=40,40,rect,yellow,yellow,2");
}

} // namespace inet

