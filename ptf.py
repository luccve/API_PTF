import math
from typing import Optional


class PTF:

    @staticmethod
    def vandenBerg(SIL: float, ARG: float, MO: float) -> float:
        carbonoOrganico = MO * 0.58
        CC_10 = (10.88 + 0.347 * ARG + 0.211 * SIL +
                 0.1756 * carbonoOrganico) / 100
        PMP = (3.83 + 0.272 * ARG + 0.212 * SIL) / 100
        return CC_10 * PMP * 10

    @staticmethod
    def arruda(SIL: float, ARG: float) -> float:
        cc_33 = (3.07439 + 0.629239 * (SIL + ARG) -
                 0.00343813 * (SIL + ARG) ** 2) / 100
        pmp = (398.889 * (SIL + ARG)) / (1308.09 + (SIL + ARG)) / 100
        return (cc_33 - pmp) * 10

    @staticmethod
    def leonorAssad(AT: float, DS: float) -> float:
        AD = (12.76278562 - 9.8726e-6 * AT ** 3) / 100
        return AD * 10

    @staticmethod
    def wenceslau(AT: float, SIL: float, ARG: float) -> float:
        ATr = AT
        SILr = SIL
        ARGr = ARG

        equaction = (
            -0.02128887 * ATr
            - 0.01005814 * SILr
            - 0.01901894 * ARGr
            + 0.0001171219 * ATr * SILr
            + 0.0002073924 * ATr * ARGr
            + 0.00006118707 * SILr * ARGr
            - 0.000006373789 * ATr * SILr * ARGr
        )

        equaction2 = equaction * 0.3591 + 1

        AD = equaction2 ** 2.78474
        return AD * 10

    @staticmethod
    def oliveira(AT: float, SIL: float, ARG: float, DS: float or None) -> float:
        densidade = DS

        if densidade == 0 or densidade is None:
            densidade = PTF.densidade(AT, ARG, None)

        AD = (-0.000021 * AT + 0.000203 * SIL + 0.000054 *
              ARG + 0.021656 * densidade) * densidade * 10
        return AD

    @staticmethod
    def barros_simplificada(AT: float, SIL: float, ARG: float) -> float:
        theta_s = 0.434714 - 0.114177 * (AT / 100) + 0.117845 * (SIL / 100)
        theta_r = 0.128617 - 0.14836 * (AT / 100) + 0.35705 * (ARG / 100)
        log_alpha = math.exp(-1.07329 - 1.59578 * (ARG / 100))
        n = 1.134153 + 0.722216 * (AT / 100) + 0.39574 * (SIL / 100)
        AD = PTF.van_genuchten(log_alpha, n, theta_r, theta_s)

        return AD

    @staticmethod
    def van_genuchten(
        log_alpha: float,
        n: float,
        theta_r: float,
        theta_s: float
    ) -> float:
        AD10 = theta_r + (theta_s - theta_r) / \
            (1 + abs(log_alpha * 10) ** n) ** (1 - 1 / n)
        AD33 = theta_r + (theta_s - theta_r) / \
            (1 + abs(log_alpha * 33) ** n) ** (1 - 1 / n)
        AD1500 = theta_r + (theta_s - theta_r) / \
            (1 + abs(log_alpha * 1500) ** n) ** (1 - 1 / n)

        return (AD33 - AD1500) * 10

    @staticmethod
    def densidade(AT: float, ARG: float, CARB: Optional[float] = None) -> float:
        bd = 1.286 + 3.208 * 0.001 * AT - 2.013 * 0.001 * ARG
        bd_oc = 1.358 + 2.79 * 0.001 * AT - 2.328 * 0.001 * ARG - 0.0052 * CARB

        if CARB is None:
            return bd
        else:
            return bd_oc

    @staticmethod
    def tomasella_n2(
        AG: float,
        AF: float,
        SIL: float,
        ARG: float,
        CARB: float,
        DS: float
    ) -> float:
        return 0.0
