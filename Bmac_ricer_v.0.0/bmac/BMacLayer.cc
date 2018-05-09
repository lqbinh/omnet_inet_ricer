//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see http://www.gnu.org/licenses/.
//

#include "inet/common/INETUtils.h"
#include "inet/common/INETMath.h"
#include "inet/networklayer/common/InterfaceEntry.h"
#include "inet/common/ModuleAccess.h"
#include "inet/linklayer/contract/IMACProtocolControlInfo.h"
#include "inet/linklayer/common/SimpleLinkLayerControlInfo.h"
#include "inet/linklayer/bmac/BMacLayer.h"

namespace inet {

Define_Module(BMacLayer);

void BMacLayer::initialize(int stage) {
    MACProtocolBase::initialize(stage);
    if (stage == INITSTAGE_LOCAL) {
        queueLength = par("queueLength");
        animation = par("animation");
        slotDuration = par("slotDuration");
        bitrate = par("bitrate");
        headerLength = par("headerLength");
        checkInterval = par("checkInterval");
        useMacAcks = par("useMACAcks");
        maxTxAttempts = par("maxTxAttempts");
        EV_DETAIL << "headerLength: " << headerLength << ", bitrate: "
                         << bitrate << endl;

        nbTxDataPackets = 0;
        nbTxPreambles = 0;
        nbRxDataPackets = 0;
        nbRxPreambles = 0;
        nbMissedAcks = 0;
        nbRecvdAcks = 0;
        nbDroppedDataPackets = 0;
        nbTxAcks = 0;

        txAttempts = 0;
        //lastDataPktDestAddr = MACAddress::BROADCAST_ADDRESS;
        //lastDataPktSrcAddr = MACAddress::BROADCAST_ADDRESS;

        macState = INIT;
        //WATCH(macState);
        nodeID = static_cast<int>(this->getAncestorPar("nodeID"));
        numSentBeacon = 0;
        numSentData = 0;

        initializeMACAddress();
        registerInterface();

        cModule *radioModule = getModuleFromPar<cModule>(par("radioModule"),
                this);
        radioModule->subscribe(IRadio::radioModeChangedSignal, this);
        radioModule->subscribe(IRadio::transmissionStateChangedSignal, this);
        radio = check_and_cast<IRadio *>(radioModule);

        wkTime = static_cast<double>(this->getAncestorPar("wkTime"));
        // init the dropped packet info
        WATCH(macState);
        WATCH(numSentBeacon);
        WATCH(numSentData);
    } else if (stage == INITSTAGE_LINK_LAYER) {

        start_ricer = new cMessage("start_ricer");
        start_ricer->setKind(START_RICER);

        wakeup = new cMessage("wakeup");
        wakeup->setKind(WAKE_UP);

        ricer_timeout = new cMessage("ricer_timeout");
        ricer_timeout->setKind(RICER_TIMEOUT);

        ricer_tx_over = new cMessage("ricer_tx_over");
        ricer_tx_over->setKind(RICER_Tx_OVER);

        ricer_stop_beacon = new cMessage("ricer_stop_beacon");
        ricer_stop_beacon->setKind(RICER_Tx_STOP_BEACON);

        ricer_stop_cca = new cMessage("ricer_stop_cca");
        ricer_stop_cca->setKind(RICER_Tx_STOP_CCA);

        ricer_stop_cbt = new cMessage("ricer_stop_cbt");
        ricer_stop_cbt->setKind(RICER_Tx_STOP_CBT);

        scheduleAt(0.0, start_ricer);
    }
}

BMacLayer::~BMacLayer() {
    cancelAndDelete(wakeup);
    cancelAndDelete(start_ricer);

    for (auto & elem : macQueue) {
        delete (elem);
    }
    macQueue.clear();
}

void BMacLayer::finish() {
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

void BMacLayer::initializeMACAddress() {
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

InterfaceEntry *BMacLayer::createInterfaceEntry() {
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
void BMacLayer::handleUpperPacket(cPacket *msg) {
    /*
     bool pktAdded = addToQueue(msg);
     if (!pktAdded)
     return;
     // force wakeup now
     if (wakeup->isScheduled() && (macState == SLEEP)) {
     cancelEvent(wakeup);
     scheduleAt(simTime() + dblrand() * 0.1f, wakeup);
     }
     */
}

void BMacLayer::sendMacBeacon() {
    BMacFrame *beacon = new BMacFrame();
    beacon->setSrcAddr(address);
    EV << "src beacon host 0 " << address << "-----------" << endl;

    beacon->setDestAddr(MACAddress::BROADCAST_ADDRESS);
    beacon->setKind(RICER_BEACON);
    beacon->setBitLength(headerLength);

    //attach signal and send down
    attachSignal(beacon);
    sendDown(beacon);
    nbTxAcks++;
}

void BMacLayer::sendMacAck(MACAddress desAddress) {
    BMacFrame *ack = new BMacFrame();
    ack->setSrcAddr(address);
    EV << "src mac Addrs host  " << address << "-----------" << endl;
    ack->setDestAddr(desAddress);
    ack->setKind(RICER_ACK);
    ack->setBitLength(headerLength);

    //attach signal and send down
    attachSignal(ack);
    sendDown(ack);
    nbTxAcks++;
}

void BMacLayer::sendMacData(MACAddress desAddress) {
    BMacFrame *data = new BMacFrame();
    EV << "src of data packet from " << address << "-----------" << endl;
    EV << "des of data packet  " << desAddress << "-----------" << endl;
    data->setSrcAddr(address);
    data->setDestAddr(desAddress);
    data->setKind(RICER_DATA);
    data->setBitLength(headerLength + 128);

    //attach signal and send down
    attachSignal(data);
    sendDown(data);
    nbTxAcks++;
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
void BMacLayer::handleSelfMessage(cMessage *msg) {

    switch (macState) {
    case INIT:
        if (msg->getKind() == START_RICER) {

            if (nodeID == 0) {
                EV_DETAIL << "This is the Base station node, beacon is sending "
                                 << endl;

                scheduleAt(simTime() + rx_wakeup_period, wakeup);

                radio->setRadioMode(IRadio::RADIO_MODE_TRANSMITTER);
                sendMacBeacon();
                numSentBeacon++;
                macState = RICER_Rx_FINISH_BEACON;
            } else {
                EV_DETAIL << "This is the Tx node, waiting for beacon " << endl;
                EV_DETAIL << "wake up time Tx "<< wkTime << endl;
                scheduleAt(simTime() + wkTime, wakeup);

                radio->setRadioMode(IRadio::RADIO_MODE_RECEIVER);
                scheduleAt(simTime() + 0.05, ricer_timeout);
                macState = RICER_Tx_IDLE_BEACON;
            }
            return;

        }
        break;

    case RICER_Tx_IDLE_BEACON:
        if (msg->getKind() == RICER_TIMEOUT) {
            EV_DETAIL << "There is no beacon, go to sleep mode " << endl;

            cancelEvent(ricer_timeout);
            radio->setRadioMode(IRadio::RADIO_MODE_SLEEP);

            macState = RICER_Tx_SLEEP;
            return;
        }
        if (msg->getKind() == RICER_BEACON) {
            EV_DETAIL << "There is a beacon, data is sending " << endl;
            cancelEvent(ricer_timeout);

//                radio->setRadioMode(IRadio::RADIO_MODE_TRANSMITTER);
//                sendMacData();
            BMacFrame *beacon = omnetpp::check_and_cast<BMacFrame *>(msg);
            lastDataPktDestAddr = beacon->getSrcAddr();

            scheduleAt(simTime(), ricer_stop_beacon);
            macState = RICER_Tx_WUB_RECV;
            return;
        }

        break;
    case RICER_Tx_WUB_RECV:
        if (msg->getKind() == RICER_Tx_STOP_BEACON) {
            EV_DETAIL << "****Tx: Beacon received... " << endl;
            //e_C += E_WAB; //WAB energy consumption
            macState = RICER_Tx_CCA;
            scheduleAt(simTime(), ricer_stop_cca);
            return;
        }
        break;
    case RICER_Tx_CCA:
        if (msg->getKind() == RICER_Tx_STOP_CCA) {
            EV_DETAIL << "****Tx: CCA done.. " << endl;
            //e_C += E_CCA; //CCA energy consumption
            macState = RICER_Tx_CBT;
            scheduleAt(simTime(), ricer_stop_cbt);
            return;
        }
        break;
    case RICER_Tx_CBT:
        if (msg->getKind() == RICER_Tx_STOP_CBT) {
            EV_DETAIL << "****Tx: CBT done.. " << endl;
            //e_C += E_CBT; //CBT idle energy consumption
            macState = RICER_Tx_FINISH_DATA;


            radio->setRadioMode(IRadio::RADIO_MODE_TRANSMITTER);
            EV << lastDataPktDestAddr
                      << " ----------------------------xoxoxo-------------------------------"
                      << endl;
            sendMacData(lastDataPktDestAddr);
            numSentData++;
            EV_DETAIL << "****Tx: Sending Mac data packet.. " << endl;
            return;
        }
        break;
    case RICER_Tx_FINISH_DATA:
        if (msg->getKind() == RICER_Tx_OVER) {
            EV_DETAIL << "Data is already sent!!!" << endl;

            scheduleAt(simTime() + 0.05, ricer_timeout);
            radio->setRadioMode(IRadio::RADIO_MODE_RECEIVER);
            macState = RICER_Tx_WAIT_ACK;

            return;
        }

        break;
    case RICER_Tx_WAIT_ACK:
        if (msg->getKind() == RICER_TIMEOUT) {
            EV_DETAIL << "There is no ack, go to sleep mode " << endl;

            cancelEvent(ricer_timeout);
            radio->setRadioMode(IRadio::RADIO_MODE_SLEEP);

            macState = RICER_Tx_SLEEP;
            return;
        }
        if (msg->getKind() == RICER_ACK) {
            EV_DETAIL << "There is an ACK" << endl;

            cancelEvent(ricer_timeout);
            radio->setRadioMode(IRadio::RADIO_MODE_SLEEP);

            macState = RICER_Tx_SLEEP;
            return;
        }
        break;
    case RICER_Tx_SLEEP:
        if (msg->getKind() == WAKE_UP) {
            EV_DETAIL << "Tx node is woken up!!!" << endl;
            scheduleAt(simTime() + wkTime, wakeup);

            radio->setRadioMode(IRadio::RADIO_MODE_RECEIVER);
            scheduleAt(simTime() + 0.05, ricer_timeout);
            macState = RICER_Tx_IDLE_BEACON;
            return;
        }
        break;
    case RICER_Rx_FINISH_BEACON:
        if (msg->getKind() == RICER_Tx_OVER) {
            EV_DETAIL << "Sending beacon finished!!!! Waiting a data..."
                             << endl;

            radio->setRadioMode(IRadio::RADIO_MODE_RECEIVER);

            macState = RICER_Rx_WAIT_DATA;
            scheduleAt(simTime() + 0.05, ricer_timeout);
            return;
        }
        break;
    case RICER_Rx_WAIT_DATA:
        if (msg->getKind() == RICER_TIMEOUT) {
            cancelEvent(ricer_timeout);
            EV_DETAIL << "There is no data, go to sleep mode" << endl;
            radio->setRadioMode(IRadio::RADIO_MODE_SLEEP);
            macState = RICER_Rx_SLEEP;

            return;
        }
        if (msg->getKind() == RICER_DATA) {
            cancelEvent(ricer_timeout);
            EV_DETAIL << "There is a data packet" << endl;

            BMacFrame *mac = static_cast<BMacFrame *>(msg);
            const MACAddress src = mac->getSrcAddr();

            EV_DETAIL << "Data packet is from : " << src << endl;

            radio->setRadioMode(IRadio::RADIO_MODE_TRANSMITTER);
            macState = RICER_Rx_FINSH_ACK;
            sendMacAck(src);

            return;
        }
        break;
    case RICER_Rx_FINSH_ACK:
        if (msg->getKind() == RICER_Tx_OVER) {
            EV_DETAIL << "ACK is sent!!!" << endl;
            radio->setRadioMode(IRadio::RADIO_MODE_SLEEP);
            macState = RICER_Rx_SLEEP;
            return;
        }
        break;
    case RICER_Rx_SLEEP:
        if (msg->getKind() == WAKE_UP) {
            EV_DETAIL << "Wake up and send beacon" << endl;
            radio->setRadioMode(IRadio::RADIO_MODE_TRANSMITTER);
            sendMacBeacon();
            numSentBeacon++;
            macState = RICER_Rx_FINISH_BEACON;
            scheduleAt(simTime() + rx_wakeup_period, wakeup);
            return;
        }
        break;

    }
    throw cRuntimeError(
            "Undefined event of type %d in state %d (radio mode %d, radio reception state %d, radio transmission state %d)!",
            msg->getKind(), macState, radio->getRadioMode(),
            radio->getReceptionState(), radio->getTransmissionState());
}

/**
 * Handle BMAC preambles and received data packets.
 */
void BMacLayer::handleLowerPacket(cPacket *msg) {
    handleSelfMessage(msg);
}

void BMacLayer::receiveSignal(cComponent *source, simsignal_t signalID,
        long value, cObject *details) {
    Enter_Method_Silent();

    // Transmission of one packet is over
    if (signalID == IRadio::transmissionStateChangedSignal) {
        IRadio::TransmissionState newRadioTransmissionState = (IRadio::TransmissionState)value;
        if (transmissionState == IRadio::TRANSMISSION_STATE_TRANSMITTING && newRadioTransmissionState == IRadio::TRANSMISSION_STATE_IDLE) {
            //PHY send a signal: Tx OVER
            scheduleAt(simTime(), ricer_tx_over);
        }
        transmissionState = newRadioTransmissionState;
    }

}

/**
 * Encapsulates the received network-layer packet into a BMacFrame and set all
 * needed header fields.
 */
bool BMacLayer::addToQueue(cMessage *msg) {
    if (macQueue.size() >= queueLength) {
        // queue is full, message has to be deleted
        EV_DETAIL << "New packet arrived, but queue is FULL, so new packet is"
                " deleted\n";
        emit(packetFromUpperDroppedSignal, msg);
        nbDroppedDataPackets++;
        return false;
    }

    BMacFrame *macPkt = encapsMsg((cPacket *) msg);
    macQueue.push_back(macPkt);
    EV_DETAIL << "Max queue length: " << queueLength << ", packet put in queue"
            "\n  queue size: " << macQueue.size() << " macState: " << macState
                     << endl;
    return true;
}

void BMacLayer::flushQueue() {
    // TODO:
    macQueue.clear();
}

void BMacLayer::clearQueue() {
    macQueue.clear();
}

void BMacLayer::attachSignal(BMacFrame *macPkt) {
    //calc signal duration
    simtime_t duration = macPkt->getBitLength() / bitrate;
    //create and initialize control info with new signal
    macPkt->setDuration(duration);
}

void BMacLayer::setColor(int color) {
    cDisplayString& dispStr = this->getDisplayString();
    if (color == 0) {
        dispStr.setTagArg("b", 3, "black");
    } else if (color == 1) //rx
            {
        dispStr.setTagArg("b", 3, "red");
    } else if (color == 2) //tx
            {
        dispStr.setTagArg("b", 3, "green");
    }
}

/**
 * Change the color of the node for animation purposes.
 */
void BMacLayer::refreshDisplay() const {
    if (!animation)
        return;
    cDisplayString& dispStr = findContainingNode(this)->getDisplayString();

    if (nodeID == 0) {
        char buf[40];
        sprintf(buf, "number of send beacon: %ld", numSentBeacon);
        dispStr.setTagArg("t", 0, buf);
    } else {
        char buf[40];
        sprintf(buf, "number of send data: %ld", numSentData);
        dispStr.setTagArg("t", 0, buf);
    }

    if (radio->getRadioMode() == IRadio::RADIO_MODE_SLEEP) {
        dispStr.setTagArg("b", 3, "yellow");

    } else if (radio->getRadioMode() == IRadio::RADIO_MODE_RECEIVER) {
        dispStr.setTagArg("b", 3, "red");
    } else if (radio->getRadioMode() == IRadio::RADIO_MODE_TRANSMITTER) {
        dispStr.setTagArg("b", 3, "green");
    } else {
        dispStr.setTagArg("b", 3, "");
    }
}

cPacket *BMacLayer::decapsMsg(BMacFrame *msg) {
    cPacket *m = msg->decapsulate();
    setUpControlInfo(m, msg->getSrcAddr());
    // delete the macPkt
    delete msg;
    EV_DETAIL << " message decapsulated " << endl;
    return m;
}

BMacFrame *BMacLayer::encapsMsg(cPacket *netwPkt) {
    BMacFrame *pkt = new BMacFrame(netwPkt->getName(), netwPkt->getKind());
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
cObject *BMacLayer::setUpControlInfo(cMessage * const pMsg,
        const MACAddress& pSrcAddr) {
    SimpleLinkLayerControlInfo * const cCtrlInfo =
            new SimpleLinkLayerControlInfo();
    cCtrlInfo->setSrc(pSrcAddr);
    cCtrlInfo->setInterfaceId(interfaceEntry->getInterfaceId());
    pMsg->setControlInfo(cCtrlInfo);
    return cCtrlInfo;
}

} // namespace inet

