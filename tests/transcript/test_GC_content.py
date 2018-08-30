"""
Test GC_content
"""

from tests.transcript.config import transcript2


def test_gc_content_CDS(transcript2):
    res = transcript2.gc_content(0)
    assert res == 0.4835164835164835


def test_gc_content_utr5(transcript2):
    res = transcript2.gc_content(1)
    assert res == 0.5769230769230769


def test_gc_content_utr3(transcript2):
    res = transcript2.gc_content(2)
    assert res == 0.3631361760660248
