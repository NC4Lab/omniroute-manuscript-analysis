"""
Module: io_rosbag.py

Purpose:
    Utilities for extracting ROS data from .bag files.
"""

from pathlib import Path
import numpy as np
from rosbags.highlevel import AnyReader as RosBagReader

from preprocessing.utils.omni_anal_logger import omni_anal_logger


def load_ros_sync_ts(bag_path: Path) -> np.ndarray:
    """
    Extract sync pulse timestamps from the '/event' topic in a ROS2 .bag file.
    Sync pulses are identified by Event messages with label 'sync_spikegadgets'.

    Parameters:
        bag_path (Path): Path to a ROS2 .bag file.

    Returns:
        np.ndarray: Array of sync pulse times (in seconds).
    """
    
    SYNC_TOPIC = "/event"
    SYNC_LABEL = "sync_spikegadgets"
    ros_sync_ts = []

    with RosBagReader([bag_path]) as reader:
        connections = [x for x in reader.connections if x.topic == SYNC_TOPIC]
        if not connections:
            raise ValueError(f"No connections found for topic '{SYNC_TOPIC}' in {bag_path}")

        for connection, timestamp, rawdata in reader.messages(connections=connections):
            msg = reader.deserialize(rawdata, connection.msgtype)

            # Only collect timestamps from sync_spikegadgets events
            if hasattr(msg, "event") and str(msg.event) == SYNC_LABEL:
                ros_sync_ts.append(timestamp / 1e9)

    if not ros_sync_ts:
        raise ValueError(f"No sync events labeled '{SYNC_LABEL}' found in topic '{SYNC_TOPIC}' from {bag_path}")
    else:
        omni_anal_logger.info(f"Extracted {len(ros_sync_ts)} sync events from topic '{SYNC_TOPIC}' with label '{SYNC_LABEL}'")

    return np.array(ros_sync_ts)


def print_all_topics(bag_path: Path) -> None:
    """
    Print all topics available in the given ROS2 .bag file.

    Parameters:
        bag_path (Path): Path to the ROS .bag file.
    """
    from rosbags.highlevel import AnyReader as RosBagReader

    with RosBagReader([bag_path]) as reader:
        topic_list = sorted(set(conn.topic for conn in reader.connections))
        if not topic_list:
            omni_anal_logger.info(f"No topics found in: {bag_path}")
        else:
            omni_anal_logger.info(f"Topics in {bag_path}:")
            for topic in topic_list:
                omni_anal_logger.info(f"  - {topic}")


def print_unique_strings_from_topic(bag_path: Path, topic: str) -> None:
    """
    Print all unique string values found in messages on the specified topic.

    Parameters:
        bag_path (Path): Path to the ROS2 .bag file.
        topic (str): Name of the topic to inspect.
    """

    unique_strings = set()

    with RosBagReader([bag_path]) as reader:
        connections = [x for x in reader.connections if x.topic == topic]
        if not connections:
            omni_anal_logger.info(f"No connections found for topic '{topic}' in {bag_path}")
            return

        for connection, timestamp, rawdata in reader.messages(connections=connections):
            msg = reader.deserialize(rawdata, connection.msgtype)

            # Case 1: The message itself is a string (e.g., std_msgs/String)
            if isinstance(msg, str):
                unique_strings.add(msg)

            # Case 2: Message has a 'data' or 'label' field that is a string
            elif hasattr(msg, "data") and isinstance(msg.data, str):
                unique_strings.add(msg.data)
            elif hasattr(msg, "label") and isinstance(msg.label, str):
                unique_strings.add(msg.label)

            # Case 3: Fallback — look through all fields
            else:
                for attr in dir(msg):
                    if not attr.startswith("_"):
                        val = getattr(msg, attr)
                        if isinstance(val, str):
                            unique_strings.add(val)

    if not unique_strings:
        omni_anal_logger.info(f"No string values found in topic '{topic}'")
    else:
        omni_anal_logger.info(f"Unique strings found in topic '{topic}':")
        for s in sorted(unique_strings):
            omni_anal_logger.info(f"  - '{s}'")
