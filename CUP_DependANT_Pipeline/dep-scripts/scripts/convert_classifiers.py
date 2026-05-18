from vulnerability_classifier import VulnerabilityClassifier


def convert(filename):
    vc = VulnerabilityClassifier.load(filename)
    vc.save(filename)
