# Failure Log

Record of what went wrong, why, and what we learned.

## F1: Wiley book chapter downloads blocked by Cloudflare Turnstile
**Date:** 2026-03-24
**What:** wget, curl, and headless Playwright all fail to download chapter PDFs from Wiley Online Library, despite TIB VPN access confirmed (book landing pages load fine).
**Why:** Cloudflare Turnstile CAPTCHA on individual chapter/PDF URLs. Even with webdriver spoofing, the Turnstile challenge blocks automated download.
**Workaround:** Book landing pages load in stealth headless mode. PDFs must be downloaded manually via browser. 4 books affected: Rybicki & Lightman, Shapiro & Teukolsky, Mészáros, Mihalas.
**Impact:** The key equations from these books are also available in the arXiv papers we successfully downloaded (Potekhin papers contain the opacity formulae; Haakonsen describes the Feautrier method; Potekhin 2013 has the TOV/EOS). Books needed primarily for additional verification and the original derivations.
